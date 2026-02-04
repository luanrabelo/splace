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
    
    - Returns:
        **str**: Path to output file if successful, None otherwise
    """
    input_gb_file = kwargs.get("genbank_file", None)
    output_fasta_path = kwargs.get("output_path", None)
    data_type = kwargs.get("data_type", None)
    executor = kwargs.get("executor", None)

    if not input_gb_file or not os.path.exists(input_gb_file):
        logging.error(f"Input GenBank file ({input_gb_file}) not found.")
        return None

    if output_fasta_path is None:
        output_fasta_path = os.getcwd()
    
    os.makedirs(output_fasta_path, exist_ok=True)

    def _process_genbank(input_gb_file, output_fasta_path, data_type):
        sequences = []
        written_files = set()
                
        #if data_type == "mt":
        #    genes_list = ["COI", "COII", "COIII", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "ATP6", "ATP8"]
        #elif data_type == "cp":
        #    genes_list = ["rbcL", "matK", "ndhF", "atpB", "psaA", "psbA", "psbB", "psbC", "psbD", "psbE", "psbF", "psbH", "psbI", "psbJ", "psbK", "psbL", "psbM", "psbN", "psbT"]
        #else:
        #    genes_list = None  # No filtering, extract all CDS

        try:
            for record in SeqIO.parse(input_gb_file, "genbank"):
                species = getattr(record, 'annotations', {}).get('organism', f"Unknown Species in {os.path.basename(input_gb_file)}")
                
                for feature in record.features:
                    if feature.type in ("CDS"):
                        product = feature.qualifiers.get('product', [None])[0]
                        
                        if product in GENE_NAME_CACHE.values():
                            gene_name = next((k for k, v in GENE_NAME_CACHE.items() if v == product), None)
                        else:
                            gene_name = sg.fix_gene_name(geneName=product, type=str(data_type))
                        
                        if not gene_name or gene_name == "None" or gene_name == None:
                            logging.warning(msg=f"Could not standardize gene name for product '{product}' in file {os.path.basename(input_gb_file)}")
                        else:
                            if product not in GENE_NAME_CACHE.values():
                                GENE_NAME_CACHE[gene_name] = product
                            
                            output_file_path = os.path.join(output_fasta_path, f"{gene_name}.fasta")
                            with open(output_file_path, "a+") as output_file:
                                output_file.write(f">{str(species).replace(' ', '_')}_{getattr(record, 'id', 'N/A')}\n{str(feature.extract(record.seq))}\n")
                                sequences.append(feature.extract(record.seq))
                            written_files.add(output_file_path)

            
            if sequences:
                logging.info(f"Converted {len(sequences)} sequences from {os.path.basename(input_gb_file)} to FASTA")
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
            return await loop.run_in_executor(local_executor, _process_genbank, input_gb_file, output_fasta_path, data_type)
    else:
        return await loop.run_in_executor(executor, _process_genbank, input_gb_file, output_fasta_path, data_type)

async def extract_multiple_genbanks(**kwargs):
    """
    ## Convert multiple GenBank files to FASTA with concurrency control
    
    - Args:
        - `genbank_files` (**List[str]**): List of GenBank file paths
        - `output_dir` (**str**): Output directory for FASTA files
        - `mitochondrial` (**bool**): Filter for mitochondrial sequences
        - `chloroplast` (**bool**): Filter for chloroplast sequences
        - `max_concurrent` (**int**): Maximum concurrent conversions (default: 5)

    - Returns:
        **List[str]**: List of successfully converted file paths
    """
    list_gb_files = kwargs.get("genbank_files", [])
    output_directory = kwargs.get("output_dir", "markers_fasta")
    data_type = kwargs.get("data_type", None)
    max_concurrent = kwargs.get("max_concurrent", 5)
   
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
                executor=executor
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
