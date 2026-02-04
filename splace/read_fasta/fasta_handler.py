import asyncio
import logging
import os
import re

from Bio import SeqIO
from SynGenes import SynGenes
from concurrent.futures import ThreadPoolExecutor

sg = SynGenes(verbose=False)
GENE_NAME_CACHE = {}

async def fasta_extractor(**kwargs):
    """
    ## Extract sequences from FASTA file based on gene/protein attributes
    
    - Args:
        - `fasta_file` (**str**): Input FASTA file path
        - `output_path` (**str**): Output directory for extracted gene files
        - `data_type` (**str**): Type of sequences to extract ("mt" for mitochondrial, "cp" for chloroplast)
        - `executor` (**ThreadPoolExecutor**): Executor for blocking I/O
    
    - Returns:
        **List[str]**: List of paths to output files that were written to
    """
    input_fasta_file = kwargs.get("fasta_file", None)
    output_fasta_path = kwargs.get("output_path", None)
    data_type = kwargs.get("data_type", None)
    executor = kwargs.get("executor", None)

    if not input_fasta_file or not os.path.exists(input_fasta_file):
        logging.error(f"Input FASTA file ({input_fasta_file}) not found.")
        return []

    if output_fasta_path is None:
        output_fasta_path = os.getcwd()
    
    os.makedirs(output_fasta_path, exist_ok=True)

    def _process_fasta(input_file, output_path, dtype):
        sequences_written = []
        written_files = set()
        
        # Define gene lists (same as in genbank_handler)
        if dtype == "mt":
            genes_list = ["COI", "COII", "COIII", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "ATP6", "ATP8"]
        elif dtype == "cp":
            genes_list = ["rbcL", "matK", "ndhF", "atpB", "psaA", "psbA", "psbB", "psbC", "psbD", "psbE", "psbF", "psbH", "psbI", "psbJ", "psbK", "psbL", "psbM", "psbN", "psbT"]
        else:
            genes_list = None

        # Base name for the header
        file_base_name = os.path.splitext(os.path.basename(input_file))[0]

        try:
            for record in SeqIO.parse(input_file, "fasta"):
                description = record.description
                
                # Extract [gene=...] or [protein=...]
                gene_match = re.search(r'\[gene=([^\]]+)\]', description)
                protein_match = re.search(r'\[protein=([^\]]+)\]', description)
                
                gene_name = None
                
                # 1. Try gene attribute
                if gene_match:
                    raw_gene = gene_match.group(1)
                    # Check if it matches our target list (if filtering is enabled)
                    if genes_list:
                        # Simple case insensitive check or partial match?
                        # SynGenes normally handles this, but if we trust gene= tag:
                        if raw_gene.upper() in genes_list:
                             gene_name = raw_gene.upper()
                        else:
                             # Try to normalize if not exact match or check SynGenes
                             normalized = sg.fix_gene_name(geneName=raw_gene, type=str(dtype) if dtype else "mt") # default mt if None
                             if normalized and normalized in genes_list:
                                 gene_name = normalized
                    else:
                        gene_name = raw_gene # No filter, extract everything

                # 2. If no gene name yet, try checking the description header manually
                if not gene_name:
                    parts = description.split()
                    # format: >atp8_ITV1046I2 atp8 ATP synthase F0 subunit 8 7816:7956 forward
                    if len(parts) > 1:
                        candidate_gene = parts[1]
                        # Check if simple candidate is in the list
                        if genes_list and candidate_gene.upper() in genes_list:
                            gene_name = candidate_gene.upper()
                        elif genes_list:
                            # Try SynGenes on the candidate token
                            normalized = sg.fix_gene_name(geneName=candidate_gene, type=str(dtype) if dtype else "mt")
                            if normalized and normalized in genes_list:
                                gene_name = normalized
                            else:
                                # Construct full name excluding coordinates (digit:digit)
                                valid_parts = []
                                for p in parts[2:]:
                                    if re.match(r'\d+:\d+', p):
                                        break
                                    valid_parts.append(p)
                                
                                if valid_parts:
                                    full_header_protein = " ".join(valid_parts)
                                    # Try normalizing the full description string
                                    normalized_full = sg.fix_gene_name(geneName=full_header_protein, type=str(dtype) if dtype else "mt")
                                    if normalized_full and normalized_full in genes_list:
                                        gene_name = normalized_full

                # 3. If no gene name yet, try protein attribute
                if not gene_name and protein_match:
                    raw_protein = protein_match.group(1)
                    
                    if raw_protein in GENE_NAME_CACHE.values():
                        gene_name = next((k for k, v in GENE_NAME_CACHE.items() if v == raw_protein), None)
                    else:
                        normalized = sg.fix_gene_name(geneName=raw_protein, type=str(dtype) if dtype else "mt")
                        if normalized:
                            gene_name = normalized
                            if raw_protein not in GENE_NAME_CACHE.values():
                                GENE_NAME_CACHE[gene_name] = raw_protein
                
                # 3. Final check against filter list
                if gene_name:
                    if genes_list and gene_name not in genes_list:
                        # Skip if it is a gene but not one we want
                        continue
                    
                    # Write sequence
                    output_file_path = os.path.join(output_path, f"{gene_name}.fasta")
                    with open(output_file_path, "a+") as out_f:
                        # Use file name as header ID + gene info for traceability if needed? 
                        # User requested: "o cabeçãlho pode ser o nome do arquivo pra facilitar"
                        out_f.write(f">{file_base_name}\n{str(record.seq)}\n")
                    
                    written_files.add(output_file_path)
                    sequences_written.append(gene_name)

            if sequences_written:
                logging.info(f"Extracted {len(sequences_written)} sequences from {os.path.basename(input_file)}")
                return list(written_files)
            else:
                logging.warning(f"No matching genes found in {os.path.basename(input_file)}")
                return []

        except Exception as e:
            logging.error(f"Error processing FASTA {os.path.basename(input_file)}: {e}")
            return []

    loop = asyncio.get_running_loop()
    if executor is None:
        with ThreadPoolExecutor() as local_executor:
            return await loop.run_in_executor(local_executor, _process_fasta, input_fasta_file, output_fasta_path, data_type)
    else:
        return await loop.run_in_executor(executor, _process_fasta, input_fasta_file, output_fasta_path, data_type)

async def extract_multiple_fastas(**kwargs):
    """
    ## Process multiple FASTA files
    """
    list_files = kwargs.get("fasta_files", [])
    output_directory = kwargs.get("output_dir", "markers_fasta")
    data_type = kwargs.get("data_type", None)
    max_concurrent = kwargs.get("max_concurrent", 5)

    if not list_files:
        return []

    if not os.path.exists(output_directory):
        os.makedirs(output_directory, exist_ok=True, mode=0o755)

    semaphore = asyncio.Semaphore(max_concurrent)
    executor = ThreadPoolExecutor(max_workers=max_concurrent)

    async def extract_single_file(f_file):
        async with semaphore:
            return await fasta_extractor(
                fasta_file=f_file,
                output_path=output_directory,
                data_type=data_type,
                executor=executor
            )

    tasks = [extract_single_file(f) for f in list_files]
    logging.info(f"Starting extraction from {len(list_files)} FASTA files")
    
    results = await asyncio.gather(*tasks, return_exceptions=True)
    executor.shutdown(wait=False)

    extracted_files = set()
    for i, result in enumerate(results):
        if isinstance(result, Exception):
            logging.error(f"Error processing {list_files[i]}: {result}")
        elif isinstance(result, list):
            extracted_files.update(result)
    
    logging.info(f"FASTA extraction completed. {len(extracted_files)} unique gene files updated.")
    
    return list(extracted_files)
