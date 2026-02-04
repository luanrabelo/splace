import asyncio
import logging
import os
import shutil
import sys

async def trimal_trimming(**kwargs):
    """
    ## Run TrimAl on a single aligned FASTA file
    
    - Args:
        - `input_file` (**str**): Input aligned FASTA file path
        - `output_file` (**str**): Output trimmed FASTA file path
        - `trimal_params` (**str**): TrimAl parameters (default: "-automated1")
        - `timeout` (**int**): Timeout in seconds (default: 3600)
    
    - Returns:
        **str**: Path to output file if successful, None otherwise
    """
    input_file = kwargs.get("input_file", None)
    output_file = kwargs.get("output_file", None)
    trimal_params = kwargs.get("trimal_params", "-automated1")
    timeout = kwargs.get("timeout", 3600)

    if shutil.which("trimal") is None:
        logging.error("TrimAl is not installed or not found in PATH. Please install TrimAl and try again.")
        sys.exit(1)
    
    if not input_file or not os.path.exists(input_file):
        logging.error(f"Input aligned file ({input_file}) not found.")
        return None
    
    if output_file is None:
        base_name = os.path.splitext(input_file)[0]
        output_file = f"{base_name}_trimmed.fasta"
    
    cmd = ["trimal", "-in", input_file, "-out", output_file]

    if trimal_params:
        cmd.extend(trimal_params.split())

    try:
        logging.info(f"Starting TrimAl for {os.path.basename(input_file)}, please wait...")

        process = await asyncio.create_subprocess_exec(
            *cmd,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE
        )
        
        try:
            stdout, stderr = await asyncio.wait_for(
                process.communicate(),
                timeout=timeout
            )
        except asyncio.TimeoutError:
            process.kill()
            await process.wait()
            logging.error(f"TrimAl timed out for {os.path.basename(input_file)}")
            return None
        
        if process.returncode == 0:
            logging.info(f"TrimAl completed for {os.path.basename(output_file)}")
            return output_file
        else:
            error_msg = stderr.decode('utf-8')
            msg = error_msg if error_msg else stdout.decode('utf-8')
            logging.error(f"TrimAl error for {os.path.basename(input_file)}: {msg}")
            return None
            
    except Exception as e:
        logging.error(f"Error running TrimAl for {os.path.basename(input_file)}: {e}")
        return None

async def trim_multiple_files(**kwargs):
    """
    ## Trim multiple aligned FASTA files using TrimAl with concurrency control
    
    - Args:
        - `alignment_files` (**List[str]**): List of aligned FASTA file paths
        - `output_dir` (**str**): Output directory for trimmed files
        - `max_concurrent` (**int**): Maximum concurrent processes (default: 5)
        - `trimal_params` (**str**): TrimAl parameters
        - `timeout` (**int**): Timeout per process

    - Returns:
        **List[str]**: List of successfully trimmed file paths
    """
    list_files = kwargs.get("alignment_files", [])
    output_directory = kwargs.get("output_dir", "trimmed_markers")
    max_concurrent = kwargs.get("max_concurrent", 5)
    trimal_params = kwargs.get("trimal_params", "-automated1")
    timeout = kwargs.get("timeout", 3600)
    
    if not list_files:
        logging.warning("No alignment files provided for trimming")
        return []
    
    if not os.path.exists(output_directory):
        logging.info(f"Creating output directory: {output_directory}")
        os.makedirs(output_directory, exist_ok=True, mode=0o755)
    
    semaphore = asyncio.Semaphore(max_concurrent)
    
    async def trim_single_file(**kwargs):
        input_file = kwargs.get('input_file')
        async with semaphore:
            base_name = os.path.splitext(os.path.basename(input_file))[0]
            output_file = os.path.join(output_directory, f"{base_name}_trimmed.fasta")
            
            return await trimal_trimming(
                input_file=input_file,
                output_file=output_file,
                trimal_params=trimal_params,
                timeout=timeout
            )
    
    tasks = [trim_single_file(input_file=f) for f in list_files]
    logging.info(f"Starting trimming of {len(list_files)} files with {max_concurrent} concurrent processes")
    results = await asyncio.gather(*tasks, return_exceptions=True)
    
    trimmed_files = []
    for i, result in enumerate(results):
        if isinstance(result, Exception):
            logging.error(f"Error trimming {list_files[i]}: {result}")
        elif result is not None:
            trimmed_files.append(result)
    
    logging.info(f"Trimming completed: {len(trimmed_files)}/{len(list_files)} successful")
    
    return trimmed_files
