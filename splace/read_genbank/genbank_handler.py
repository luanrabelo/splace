import asyncio
import logging
import os

from Bio import SeqIO
from SynGenes import SynGenes

from concurrent.futures import ThreadPoolExecutor

GENE_NAME_CACHE = {}
sg = SynGenes(verbose=False)

async def genbank_to_fasta(**kwargs):
    """
    ## Convert GenBank file to FASTA format

    - Args:
        - `genbank_file` (**str**): Input GenBank file path
        - `output_path` (**str**): Output FASTA file path
        - `data_type` (**str**): Type of sequences to extract ("mt" for mitochondrial, "cp" for chloroplast)
        - `executor` (**ThreadPoolExecutor**): Executor for blocking I/O
        - `genes_filter` (**List[str]**, optional): List of specific gene names to extract. Overrides default genes_list.
        - `feature_types` (**List[str]**, optional): List of GenBank feature types to scan (e.g., ["CDS", "rRNA", "tRNA"]). Default: ["CDS"].

    - Returns:
        **str**: Path to output file if successful, None otherwise
    """
    input_gb_file = kwargs.get("genbank_file", None)
    output_fasta_path = kwargs.get("output_path", None)
    data_type = kwargs.get("data_type", None)
    executor = kwargs.get("executor", None)
    genes_filter = kwargs.get("genes_filter", None)
    feature_types = kwargs.get("feature_types", None)

    if not input_gb_file or not os.path.exists(input_gb_file):
        logging.error(f"Input GenBank file ({input_gb_file}) not found.")
        return None

    if output_fasta_path is None:
        output_fasta_path = os.getcwd()
    
    os.makedirs(output_fasta_path, exist_ok=True)

    def _process_genbank(input_gb_file, output_fasta_path, data_type, genes_filter, feature_types):
        written_files = set()

        # Determine which feature types to scan
        if feature_types:
            target_feature_types = feature_types
        else:
            target_feature_types = ["CDS"]

        # Determine gene name filter list
        if genes_filter:
            genes_list = genes_filter
        elif data_type == "mt":
            genes_list = ["COI", "COII", "COIII", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "ATP6", "ATP8"]
        elif data_type == "cp":
            genes_list = ["rbcL", "matK", "ndhF", "atpB", "psaA", "psbA", "psbB", "psbC", "psbD", "psbE", "psbF", "psbH", "psbI", "psbJ", "psbK", "psbL", "psbM", "psbN", "psbT"]
        else:
            genes_list = None  # No filtering, extract all

        # Collect best (longest) sequence per gene per species to handle duplicates
        best_per_gene = {}  # (gene_name, species_key) -> (header, seq_str, seq_len)

        try:
            for record in SeqIO.parse(input_gb_file, "genbank"):
                species = getattr(record, 'annotations', {}).get('organism', f"Unknown Species in {os.path.basename(input_gb_file)}")

                for feature in record.features:
                    if feature.type in target_feature_types:
                        product = feature.qualifiers.get('product', [None])[0]
                        gene_qualifier = feature.qualifiers.get('gene', [None])[0]

                        # Try gene qualifier first, then product
                        raw_name = gene_qualifier or product

                        if not raw_name:
                            continue

                        if raw_name in GENE_NAME_CACHE:
                            gene_name = GENE_NAME_CACHE[raw_name]
                        else:
                            gene_name = sg.fix_gene_name(geneName=raw_name, type=str(data_type)) if data_type else raw_name
                            if gene_name and gene_name != "None":
                                GENE_NAME_CACHE[raw_name] = gene_name

                        if not gene_name or gene_name == "None" or gene_name is None:
                            logging.warning(msg=f"Could not standardize gene name for '{raw_name}' in file {os.path.basename(input_gb_file)}")
                            continue

                        # Filter by genes_list if defined
                        if genes_list and gene_name not in genes_list:
                            continue

                        header = f">{str(species).replace(' ', '_')}_{getattr(record, 'id', 'N/A')}"
                        seq_str = str(feature.extract(record.seq))
                        seq_len = len(seq_str)
                        species_key = f"{species}_{record.id}"
                        key = (gene_name, species_key)

                        if key not in best_per_gene or seq_len > best_per_gene[key][2]:
                            if key in best_per_gene:
                                logging.info(f"Duplicate gene '{gene_name}' in {species_key}: keeping longer ({seq_len}bp over {best_per_gene[key][2]}bp)")
                            best_per_gene[key] = (header, seq_str, seq_len)

            # Write best sequences to files
            for (gene_name, _), (header, seq_str, _) in best_per_gene.items():
                output_file_path = os.path.join(output_fasta_path, f"{gene_name}.fasta")
                with open(output_file_path, "a+") as output_file:
                    output_file.write(f"{header}\n{seq_str}\n")
                written_files.add(output_file_path)

            if best_per_gene:
                logging.info(f"Converted {len(best_per_gene)} sequences from {os.path.basename(input_gb_file)} to FASTA")
                return list(written_files)
            else:
                logging.warning(f"No matching sequences found in {os.path.basename(input_gb_file)}")
                return []

        except Exception as e:
            logging.error(f"Error converting {os.path.basename(input_gb_file)}: {e}")
            return []

    loop = asyncio.get_running_loop()
    if executor is None:
        # Create a local executor if not provided
        with ThreadPoolExecutor() as local_executor:
            return await loop.run_in_executor(local_executor, _process_genbank, input_gb_file, output_fasta_path, data_type, genes_filter, feature_types)
    else:
        return await loop.run_in_executor(executor, _process_genbank, input_gb_file, output_fasta_path, data_type, genes_filter, feature_types)

async def extract_multiple_genbanks(**kwargs):
    """
    ## Convert multiple GenBank files to FASTA with concurrency control

    - Args:
        - `genbank_files` (**List[str]**): List of GenBank file paths
        - `output_dir` (**str**): Output directory for FASTA files
        - `data_type` (**str**): Type filter ("mt" or "cp")
        - `max_concurrent` (**int**): Maximum concurrent conversions (default: 5)
        - `genes_filter` (**List[str]**, optional): List of specific gene names to extract
        - `feature_types` (**List[str]**, optional): List of GenBank feature types to scan

    - Returns:
        **List[str]**: List of successfully converted file paths
    """
    list_gb_files = kwargs.get("genbank_files", [])
    output_directory = kwargs.get("output_dir", "markers_fasta")
    data_type = kwargs.get("data_type", None)
    max_concurrent = kwargs.get("max_concurrent", 5)
    genes_filter = kwargs.get("genes_filter", None)
    feature_types = kwargs.get("feature_types", None)
   
    if not list_gb_files:
        logging.warning("No GenBank files provided for conversion")
        return []

    if not os.path.exists(output_directory):
        logging.info(f"Creating output directory: {output_directory}")
        os.makedirs(output_directory, exist_ok=True, mode=0o755)

    semaphore = asyncio.Semaphore(max_concurrent)
    executor = ThreadPoolExecutor(max_workers=max_concurrent)

    async def convert_single_file(gb_file):
        async with semaphore:
            return await genbank_to_fasta(
                genbank_file=gb_file,
                output_path=output_directory,
                data_type=data_type,
                executor=executor,
                genes_filter=genes_filter,
                feature_types=feature_types
            )

    tasks = [convert_single_file(gb_file) for gb_file in list_gb_files]
    logging.info(f"Starting conversion of {len(list_gb_files)} GenBank files")
    
    results = await asyncio.gather(*tasks, return_exceptions=True)
    executor.shutdown(wait=False)

    converted_files = set()
    for i, result in enumerate(results):
        if isinstance(result, Exception):
            logging.error(f"Error converting {list_gb_files[i]}: {result}")
        elif isinstance(result, list):
            converted_files.update(result)
        elif result:
             # Fallback for unexpected single returns
             converted_files.add(result)
            
    logging.info(f"Conversion completed. {len(converted_files)} unique gene files generated/updated.")
    
    return list(converted_files)
