#!/usr/bin/env python3

import asyncio
import argparse
import logging
import shutil
import subprocess
import os
import sys

__authors__     = "Renato Oliveira and Luan Rabelo"
__license__     = "GPL-3.0"
__version__     = "4.0.0"
__maintainers__ = "Renato Oliveira and Luan Rabelo"
__email__       = "luan.rabelo@pq.itv.org"
__date__        = "2025/11/15"
__github__      = "itvgenomics/splace"
__status__      = "Stable"
__tool__        = "SPLACE"


logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y/%m/%d - %H:%M:%S'
    )

parser = argparse.ArgumentParser(
    prog=__tool__,
    usage='%(prog)s [options]',
    #description=pipenote_ascii_art,
    epilog=f"In Development...",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    add_help=True,
    )

input_group = parser.add_argument_group(
    title="Input Files",
    description="Specify the directory that contains the sequence files to be processed."
)
input_group.add_argument(
    "-i", "--input_dir",
    required=True,
    type=str,
    help="Directory with FASTA (.fasta, .fa, .fas) or GenBank (.gb, .gbk) files. Default: 'input_files/'.",
    default="input_files/"
)

output_group = parser.add_argument_group(
    title="Output Directory",
    description="Specify the directory where output files will be saved."
)
output_group.add_argument(
    "-o", "--output_dir",
    required=False,
    type=str,
    help="Directory to save output files. Default: 'splace_output/'.",
    default="splace_output/"
)

threads = parser.add_argument_group(
    title="Threads",
    description="Number of threads to use."
    )
threads.add_argument(
    "-t", "--threads",
    default=8,
    type=int,
    help="Number of threads to use. Default is 8."
    )

genbank_group = parser.add_argument_group(
    title="GenBank Options",
    description="Options for processing GenBank files."
)
genbank_group.add_argument(
    "--gb-type",
    choices=["mt", "cp"],
    help="Extract mitochondrial or chloroplast sequences from GenBank files.",
    default="mt"
)

pipeline_group = parser.add_argument_group(
    title="Pipeline Options",
    description="Options for controlling the analysis pipeline."
)
pipeline_group.add_argument(
    "--align",
    action="store_true",
    help="Perform multiple sequence alignment using MAFFT."
)
pipeline_group.add_argument(
    "--trimal",
    action="store_true",
    help="Perform alignment trimming using TrimAl. Requires --align."
)
pipeline_group.add_argument(
    "--iqtree",
    action="store_true",
    help="Perform phylogenetic analysis using IQ-TREE. Requires --trimal."
)
pipeline_group.add_argument(
    "--benchmark",
    action="store_true",
    help="Enable benchmarking of execution time."
)

if __name__ == "__main__":
    print(f"\n{'#'*70}\n")
    print(f"Version:     {__version__}")
    print(f"Status:      {__status__}")
    print(f"Authors:     {__authors__}")
    print(f"Maintainers: {__maintainers__}")
    print(f"Contact:     {__email__}")
    print(f"License:     {__license__}")
    print(f"GitHub:      {__github__}")
    print(f"\n{'#'*70}\n")
    
    tool_env = "splace_env"
    conda_env = os.getenv("CONDA_DEFAULT_ENV")

    if not shutil.which("conda"):
        logging.error("Conda is not installed or not found in PATH. Please install Conda and try again.")
        sys.exit(1)

    try:
        conda_envs = subprocess.run(
            ["conda", "env", "list"],
            capture_output=True,
            text=True,
            check=True
        ).stdout.splitlines()
    except subprocess.CalledProcessError:
        logging.error("Failed to list Conda environments.")
        sys.exit(1)
    
    if not any(tool_env in env for env in conda_envs):
        logging.error(f"The '{tool_env}' conda environment does not exist. Please create it using the provided environment.yml file.")
        sys.exit(1)

    if conda_env != tool_env:
        logging.warning(f"Activate the '{tool_env}' conda environment before running the {__tool__}. Current environment: {conda_env}. Run 'conda activate {tool_env}' and try again.")
        sys.exit(1)

    args = parser.parse_args()

    # Initialize Benchmark
    from splace.utils import Benchmark
    benchmark = Benchmark(output_path=os.path.join(args.output_dir, "benchmark.tsv"), enabled=args.benchmark)

    if not os.path.exists(args.input_dir):
        logging.error(msg=f"Input directory '{args.input_dir}' not found.")
        sys.exit(1)
        
    os.makedirs(name=args.output_dir, exist_ok=True, mode=0o755)

    genbank_files = []
    fasta_files = []
    
    for root, dirs, files in os.walk(args.input_dir):
        for file in files:
            if file.lower().endswith(('.gb', '.gbk', '.genbank')):
                genbank_files.append(os.path.join(root, file))
    
    for root, dirs, files in os.walk(args.input_dir):
        for file in files:
            if file.lower().endswith(('.fasta', '.fa', '.fas')):
                fasta_files.append(os.path.join(root, file))
    
    logging.info(msg=f"Found {len(fasta_files)} FASTA files and {len(genbank_files)} GenBank files in '{args.input_dir}'.")
    
    if genbank_files or fasta_files:
        from splace import extract_multiple_genbanks
        
        if args.trimal and not args.align:
            logging.error("The argument --trimal requires --align to be specified.")
            sys.exit(1)
        
        if args.iqtree and not args.trimal:
            logging.error("The argument --iqtree requires --trimal to be specified.")
            sys.exit(1)

        if os.path.exists(os.path.join(args.output_dir, "markers_fasta")):
            logging.error(msg=f"Output directory '{os.path.join(args.output_dir, 'markers_fasta')}' already exists. Please remove it and try again.")
            sys.exit(1)

        if os.path.exists(os.path.join(args.output_dir, "aligned_markers")):
            logging.error(msg=f"Output directory '{os.path.join(args.output_dir, 'aligned_markers')}' already exists. Please remove it and try again.")
            sys.exit(1)
        
        if os.path.exists(os.path.join(args.output_dir, "phylogeny")):
            logging.error(msg=f"Output directory '{os.path.join(args.output_dir, 'phylogeny')}' already exists. Please remove it and try again.")
            sys.exit(1)


        logging.info(msg=f"Found {len(genbank_files)} GenBank files and {len(fasta_files)} FASTA files. Processing...")
        try:
            converted_files = set()
            
            if genbank_files:
                benchmark.start("GenBank Extraction")
                gb_results = asyncio.run(
                    extract_multiple_genbanks(
                        genbank_files=genbank_files,
                        output_dir=os.path.join(args.output_dir, "markers_fasta"),
                        data_type=args.gb_type,
                        max_concurrent=args.threads
                    )
                )
                benchmark.stop("GenBank Extraction")
                if gb_results:
                    converted_files.update(gb_results)

            if fasta_files:
                from splace import extract_multiple_fastas
                benchmark.start("FASTA Extraction")
                fasta_results = asyncio.run(
                    extract_multiple_fastas(
                        fasta_files=fasta_files,
                        output_dir=os.path.join(args.output_dir, "markers_fasta"),
                        data_type=args.gb_type,
                        max_concurrent=args.threads
                    )
                )
                benchmark.stop("FASTA Extraction")
                if fasta_results:
                    converted_files.update(fasta_results)

            # Convert back to list for next steps
            converted_files = list(converted_files)

            if converted_files and args.align:
                from splace import align_multiple_files

                logging.info(msg=f"Starting alignment for {len(converted_files)} files...")
                benchmark.start("Alignment")
                aligned_files = asyncio.run(
                    align_multiple_files(
                        fasta_files=converted_files,
                        output_dir=os.path.join(args.output_dir, "aligned_markers"),
                        threads=args.threads
                    )
                )
                benchmark.stop("Alignment")

                if aligned_files and args.trimal:
                    from splace import trim_multiple_files

                    logging.info(msg=f"Starting trimming for {len(aligned_files)} files...")
                    benchmark.start("Trimming")
                    trimmed_files = asyncio.run(
                        trim_multiple_files(
                            alignment_files=aligned_files,
                            output_dir=os.path.join(args.output_dir, "trimmed_markers"),
                            trimal_params="-automated1",
                            max_concurrent=args.threads
                        )
                    )
                    benchmark.stop("Trimming")

                    if trimmed_files and args.iqtree:
                        from splace import run_phylogeny_pipeline
                        
                        logging.info(msg=f"Starting phylogeny analysis...")
                        benchmark.start("Phylogeny")
                        run_phylogeny_pipeline(
                            trimmed_files=trimmed_files,
                            output_dir=os.path.join(args.output_dir, "phylogeny"),
                            threads=args.threads
                        )
                        benchmark.stop("Phylogeny")
            
            # Save Benchmark Statistics
            total_input_files = len(genbank_files) + len(fasta_files)
            benchmark.save(input_dir=args.input_dir, file_count=total_input_files)

        except Exception as e:
            logging.error(msg=f"An error occurred during processing: {e}")
    