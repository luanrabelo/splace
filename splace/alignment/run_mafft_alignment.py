import asyncio
import logging
import os
import shutil
import sys

async def mafft_alignment(**kwargs):
    """
    ## Run MAFFT alignment on a single FASTA file
    
    - Args:
        - `fasta_file` (**str**): Input FASTA file path
        - `output_file` (**str**): Output aligned FASTA file path
        - `mafft_params` (**str**): MAFFT parameters
        - `threads` (**int**): Number of threads
        - `preserve_case` (**bool**): Preserve case (default: True)
        - `timeout` (**int**): Timeout in seconds (default: 3600)
    
    - Returns:
        **str**: Path to output file if successful, sys.exit(1) otherwise
    """
    input_fasta_file = kwargs.get("fasta_file", None)
    aligned_fasta_output = kwargs.get("output_file", None)
    mafft_params = kwargs.get("mafft_params", None)
    mafft_threads = kwargs.get("threads", 4)
    preserve_case = kwargs.get("preserve_case", True)
    mafft_timeout = kwargs.get("timeout", 3600)

    if shutil.which("mafft") is None:
        logging.error("MAFFT is not installed or not found in PATH. Please install MAFFT and try again.")
        sys.exit(1)
    
    if not input_fasta_file or not os.path.exists(input_fasta_file):
        logging.error(f"Input FASTA file ({os.path.abspath(input_fasta_file)}) not found. Please check the path and try again.")
        sys.exit(1)
    
    if os.path.isdir(input_fasta_file):
        logging.error(f"Input path ({os.path.abspath(input_fasta_file)}) is a directory. Please provide a valid FASTA file.")
        sys.exit(1)

    if aligned_fasta_output is None:
        base_name = os.path.splitext(input_fasta_file)[0]
        aligned_fasta_output = f"{base_name}_aligned.fasta"
    
    cmd = ["mafft", "--thread", str(mafft_threads)]

    if mafft_params:
        cmd.extend(mafft_params.split())
    
    if preserve_case:
        cmd.append("--preservecase")
    
    cmd.append(input_fasta_file)

    try:
        logging.info(f"Starting MAFFT alignment for {os.path.basename(input_fasta_file)}, please wait...")

        process = await asyncio.create_subprocess_exec(
            *cmd,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE
        )
        
        try:
            stdout, stderr = await asyncio.wait_for(
                process.communicate(),
                timeout=mafft_timeout
            )
        except asyncio.TimeoutError:
            process.kill()
            await process.wait()
            logging.error(f"MAFFT alignment timed out for {os.path.basename(input_fasta_file)}")
            sys.exit(1)
        
        if process.returncode == 0:
            with open(aligned_fasta_output, 'w', encoding='utf-8') as f:
                f.write(stdout.decode('utf-8'))
            logging.info(f"Alignment completed for {os.path.basename(aligned_fasta_output)}")
            
            return aligned_fasta_output
        
        else:
            error_msg = stderr.decode('utf-8')
            logging.error(f"MAFFT error for {os.path.basename(input_fasta_file)}: {error_msg}")
            sys.exit(1)
            
    except FileNotFoundError:
        logging.error("MAFFT not found. Please ensure it is installed and in your PATH")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Error running MAFFT for {os.path.basename(input_fasta_file)}: {e}")
        sys.exit(1)

async def align_multiple_files(**kwargs):
    """
    ## Align multiple FASTA files using MAFFT with concurrency control
    
    - Args:
        - `fasta_files` (**List[str]**): List of FASTA file paths
        - `output_dir` (**str**): Output directory for aligned files
        - `max_concurrent` (**int**): Maximum concurrent alignments (default: 5)
        - `mafft_params` (**str**): MAFFT parameters
        - `threads` (**int**): Number of threads per alignment
        - `preserve_case` (**bool**): Preserve case
        - `timeout` (**int**): Timeout per alignment

    - Returns:
        **List[str]**: List of successfully aligned file paths
    """
    list_fasta_files = kwargs.get("fasta_files", [])
    aligned_output_directory = kwargs.get("output_dir", "02-aligned_markers")
    max_concurrent_alignments = kwargs.get("max_concurrent", 15)
    mafft_params = kwargs.get("mafft_params", "--auto")
    mafft_threads = kwargs.get("threads", 4)
    preserve_case = kwargs.get("preserve_case", True)
    mafft_timeout = kwargs.get("timeout", 3600)
    
    if not list_fasta_files or len(list_fasta_files) == 0:
        logging.warning("No FASTA files provided for alignment")
        return []
    
    if not os.path.exists(aligned_output_directory):
        logging.info(f"Creating output directory: {aligned_output_directory}")
        os.makedirs(aligned_output_directory, exist_ok=True, mode=0o755)
    
    semaphore = asyncio.Semaphore(max_concurrent_alignments)
    
    async def align_single_file(**kwargs):
        """Align a single file with semaphore control"""
        fasta_file = kwargs.get('fasta_file')
        async with semaphore:
            base_name = os.path.splitext(os.path.basename(fasta_file))[0]
            output_file = os.path.join(aligned_output_directory, f"{base_name}_aligned.fasta")
            
            return await mafft_alignment(
                fasta_file=fasta_file,
                output_file=output_file,
                mafft_params=mafft_params,
                threads=mafft_threads,
                preserve_case=preserve_case,
                timeout=mafft_timeout
            )
    
    tasks = [align_single_file(fasta_file=fasta_file) for fasta_file in list_fasta_files]
    logging.info(f"Starting alignment of {len(list_fasta_files)} files with {max_concurrent_alignments} concurrent processes")
    results = await asyncio.gather(*tasks, return_exceptions=True)
    
    aligned_files = []
    for i, result in enumerate(results):
        if isinstance(result, Exception):
            logging.error(f"Error aligning {list_fasta_files[i]}: {result}")
        elif result is not None:
            aligned_files.append(result)
    
    logging.info(f"Alignments completed: {len(aligned_files)}/{len(list_fasta_files)} successful")
    
    return aligned_files