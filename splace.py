#!/usr/bin/env python3

import argparse
import logging
import shutil
import subprocess
import os
import sys

__authors__     = "Renato Oliveira and Luan Rabelo"
__license__     = "GPL-3.0"
__version__     = "2.0.0"
__maintainers__ = "Luan Rabelo"
__email__       = "luan.rabelo@pq.itv.org"
__date__        = "2025/11/15"
__github__      = "luanrabelo/splace"
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

env_type = parser.add_argument_group(
    title="Environment Type",
    description="Specify the type of environment to be used.",
    )
env_type.add_argument(
    "-et", "--env_type",
    default="conda",
    choices=["docker", "singularity", "conda"],
    type=str,
    help="Type of environment to be used. Choose between 'docker', 'singularity', or 'conda'. Default is 'conda'."
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



if __name__ == "__main__":
    print(f"\n{'#'*70}\n")
    print(f"Version:     {__version__}")
    print(f"Status:      {__status__}")
    print(f"Maintainers: {__maintainers__}")
    print(f"License:     {__license__}")
    print(f"GitHub:      {__github__}")
    print(f"\n{'#'*70}\n")
    
    tool_env = "splace_env"
    env_type = parser.parse_args().env_type

    if env_type == "conda":
        if not shutil.which("conda"):
            logging.error("Conda is not installed or not found in PATH. Please install Conda and try again.")
            sys.exit(1)

        conda_envs = subprocess.run(
            ["conda", "env", "list"],
            capture_output=True,
            text=True
            ).stdout.splitlines()
        
        if not any(tool_env in env for env in conda_envs):
            logging.error(f"The '{tool_env}' conda environment does not exist. Please create it using the provided environment.yml file.")
            sys.exit(1)

        conda_env = os.getenv('CONDA_DEFAULT_ENV')
        if conda_env != tool_env:
            logging.warning(f"Activate the '{tool_env}' conda environment before running the {__tool__}. Current environment: {conda_env}. Run 'conda activate {tool_env}' and try again.")
            sys.exit(1)
    
    elif env_type == "docker":
        if not shutil.which("docker"):
            logging.error("Docker is not installed or not found in PATH. Please install Docker and try again.")
            sys.exit(1)

        logging.info("Docker is installed and found in PATH.")
        docker_images = subprocess.run(
            ["docker", "images"],
            capture_output=True,
            text=True
            ).stdout.splitlines()
        
        if not any(tool_env in image for image in docker_images):
            logging.error(f"The Docker image '{tool_env}' does not exist. Building the Docker image, please wait...")
            build_process = subprocess.run(
                ["docker", "build", "-t", tool_env, "Dockerfile/Dockerfile"],
                capture_output=True,
                text=True
            )
            
            if build_process.returncode != 0:
                logging.error("Failed to build the Docker image. Please check the Dockerfile and try again.")
                sys.exit(1)
            else:
                logging.info(f"Docker image '{tool_env}' built successfully. Testing the Docker image...")
                
                docker_test = subprocess.run(
                    ["docker", "run", "--rm", tool_env, "echo", "Docker image is working correctly."],
                    capture_output=True,
                    text=True
                )
                
                if docker_test.returncode != 0:
                    logging.error("Docker image test failed. Please check the Docker image and try again.")
                    sys.exit(1)
                else:
                    logging.info("Docker image is working correctly.")