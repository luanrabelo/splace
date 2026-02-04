import logging
import os
import shutil
import subprocess
from Bio import SeqIO

def check_iqtree():
    """Check if IQ-TREE is installed (iqtree or iqtree2)."""
    if shutil.which("iqtree2"):
        return "iqtree2"
    elif shutil.which("iqtree"):
        return "iqtree"
    else:
        return None

def prepare_supermatrix(trimmed_files, output_dir):
    """
    Concatenate trimmed alignments into a supermatrix (NEXUS format) with partition block.
    
    Args:
        trimmed_files (list): List of paths to trimmed FASTA files.
        output_dir (str): Directory to save the supermatrix.
        
    Returns:
        str: Path to the generated supermatrix file.
    """
    logging.info("Preparing supermatrix...")
    
    # Sort files to ensure deterministic order
    trimmed_files.sort()
    
    all_taxa = set()
    gene_data = [] # List of tuples: (gene_name, length, sequence_dict)
    
    # 1. Read all files
    for fasta_path in trimmed_files:
        gene_name = os.path.basename(fasta_path).replace("_aligned_trimmed.fasta", "").replace("_trimmed.fasta", "").replace(".fasta", "")
        
        records = {}
        length = 0
        try:
            for record in SeqIO.parse(fasta_path, "fasta"):
                # Use description as ID (BioPython parses >ID Description)
                # Our headers are ">Species | Accession" which BioPython might split if spaces are present
                # genbank_handler writes: f">{species} | {id}"
                # If species has spaces, SeqIO.parse might treat record.id as the first word and description as the rest.
                # It's safer to use record.description which usually contains the full header line minus the '>'
                
                # Check consistency of length
                if length == 0:
                    length = len(record.seq)
                elif len(record.seq) != length:
                    logging.warning(f"Sequence length mismatch in {gene_name}. Expected {length}, got {len(record.seq)} for {record.description}")
                
                records[record.description] = str(record.seq)
                all_taxa.add(record.description)
                
            if length > 0:
                gene_data.append((gene_name, length, records))
            else:
                logging.warning(f"Gene {gene_name} is empty or invalid. Skipping.")
                
        except Exception as e:
            logging.error(f"Error reading {fasta_path}: {e}")

    if not gene_data:
        logging.error("No valid data found to create supermatrix.")
        return None

    # 2. Build Supermatrix
    output_nexus = os.path.join(output_dir, "supermatrix.nex")
    total_chars = sum(g[1] for g in gene_data)
    num_taxa = len(all_taxa)
    
    try:
        with open(output_nexus, "w", encoding="utf-8") as out:
            out.write("#NEXUS\n")
            out.write("BEGIN DATA;\n")
            out.write(f"DIMENSIONS NTAX={num_taxa} NCHAR={total_chars};\n")
            out.write("FORMAT DATATYPE=DNA GAP=- MISSING=?;\n")
            out.write("MATRIX\n")
            
            # Sort taxa for clean output
            sorted_taxa = sorted(list(all_taxa))
            
            # Determine max name length for padding
            max_name_len = max(len(t) for t in sorted_taxa) + 2
            
            for taxon in sorted_taxa:
                # Sanitize taxon name for NEXUS (replace spaces with underscores or wrap in quotes)
                # For safety, we wrap in single quotes if it contains spaces or weird chars
                safe_taxon_name = f"'{taxon}'"
                out.write(f"{safe_taxon_name.ljust(max_name_len)} ")
                
                for _, length, records in gene_data:
                    seq = records.get(taxon, "-" * length)
                    out.write(seq)
                out.write("\n")
            
            out.write(";\nEND;\n\n")
            
            # 3. Write Partition Block (Sets)
            out.write("BEGIN SETS;\n")
            current_pos = 1
            for gene_name, length, _ in gene_data:
                end_pos = current_pos + length - 1
                out.write(f"    CHARSET {gene_name} = {current_pos}-{end_pos};\n")
                current_pos = end_pos + 1
            out.write("END;\n")
            
        logging.info(f"Supermatrix created at {output_nexus} (Taxa: {num_taxa}, Sites: {total_chars})")
        return output_nexus
        
    except Exception as e:
        logging.error(f"Failed to create supermatrix: {e}")
        return None

def run_iqtree_analysis(supermatrix_path, output_dir, threads=1):
    """
    Run IQ-TREE analysis on the supermatrix.
    
    Args:
        supermatrix_path (str): Path to the supermatrix file.
        output_dir (str): Directory for output files.
        threads (int): Number of threads.
        
    Returns:
        str: Path to the main tree file (.treefile) or None.
    """
    iqtree_cmd = check_iqtree()
    if not iqtree_cmd:
        logging.error("IQ-TREE not found in PATH. Please install IQ-TREE (iqtree or iqtree2).")
        return None # Critical failure for this step

    prefix = os.path.join(output_dir, "splace_tree")
    
    # Command: iqtree2 -s file.nex -nt T -bb 1000 -m MFP --prefix ...
    # -s: Input (Nexus file determines alignment + partitions)
    # -nt: AUTO or number
    # -B: Ultrafast Bootstrap (1000)
    # -m: ModelFinder (MFP)
    
    cmd = [
        iqtree_cmd,
        "-s", supermatrix_path,
        "-nt", str(threads),
        "-B", "1000",
        "-m", "MFP",
        "-pre", prefix,
        "-redo" # Overwrite existing
    ]
    
    logging.info(f"Starting IQ-TREE analysis with {threads} threads...")
    logging.info(f"Command: {' '.join(cmd)}")
    
    try:
        # IQ-TREE writes to stdout/stderr. Capturing and logging might be noisy, 
        # but useful.
        process = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        if process.returncode == 0:
            logging.info("IQ-TREE analysis completed successfully.")
            tree_file = f"{prefix}.treefile"
            if os.path.exists(tree_file):
                return tree_file
            else:
                logging.warning(f"IQ-TREE finished 0 but {tree_file} not found.")
                return None
        else:
            logging.error(f"IQ-TREE failed with return code {process.returncode}")
            # Log the tail of stderr/stdout
            logging.error(f"IQ-TREE stderr: {process.stderr[-1000:]}")
            return None
            
    except Exception as e:
        logging.error(f"Error executing IQ-TREE: {e}")
        return None

def run_phylogeny_pipeline(trimmed_files, output_dir, threads=1):
    """
    Orchestrate the phylogeny step: Supermatrix -> IQ-TREE
    """
    if not trimmed_files:
        logging.warning("No input files for phylogeny.")
        return None
        
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True, mode=0o755)
        
    # Step 1: Supermatrix
    supermatrix = prepare_supermatrix(trimmed_files, output_dir)
    if not supermatrix:
        return None
        
    # Step 2: IQ-TREE
    tree_file = run_iqtree_analysis(supermatrix, output_dir, threads)
    
    return tree_file
