<p align="center">
  <img src="docs/assets/SPLACE.png" alt="SPLACE Logo" width="75%">
</p>

<p align="center">
    <h1 align="center">SPLACE (<b>SP</b>Lit, <b>A</b>lign and <b>C</b>oncatenat<b>E</b>)</h1>
</p>

# Contents Overview
- [System Overview](#system-overview)
- [Licence](#licence)
- [Getting Started](#getting-started)
  - [Prerequisites](#prerequisites)
  - [Installation](#installation)
  - [Usage](#usage)

***
&nbsp;
## System Overview
##### [:rocket: Go to Contents Overview](#contents-overview)
**SPLACE** is a Python toolkit for splitting, aligning, and concatenating gene sequences in phylogenetic pipelines. 
It integrates with **SynGenes** to standardize gene nomenclature across datasets, and uses asynchronous I/O and parallel alignment to accelerate large-scale workflows.
&nbsp;
> [!NOTE]
> This project is a fork of the original SPLACE repository, with enhancements for better performance and usability.
> See the original repository at [https://github.com/reinator/splace/](https://github.com/reinator/splace/)
&nbsp;
## Licence
##### [:rocket: Go to Contents Overview](#contents-overview)
**SPLACE** is released under the **GPL-3.0 License**. This license permits free use, modification, and distribution of the software, provided that any derivative works also adhere to the same license terms.
For more details, please see the [GPL-3.0 License](LICENSE).
&nbsp;
## Getting Started
##### [:rocket: Go to Contents Overview](#contents-overview)
### Prerequisites
Before you run **SPLACE**, make sure you have the following prerequisites installed on your system:
- **Python Environment and Package Manager**
    - Python **version 3.12 or higher**
    - conda
    - git
    - singularity or docker (optional, but recommended for containerized execution)
- **Required Software and Libraries**
    - `biopython` # For biological sequence handling and parsing
    - `syngenes`  # For gene nomenclature standardization
    - `mafft`     # For multiple sequence alignment
    - `trimal`    # For automated alignment trimming
    - `openpyxl`  # For Excel file handling
&nbsp;
### Installation
##### [:rocket: Go to Contents Overview](#contents-overview)
#### Conda
1. Download and install Miniconda or Anaconda from [https://docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html) or [https://www.anaconda.com/products/distribution](https://www.anaconda.com/products/distribution) respectively.
2. Clone the **SPLACE** repository and install the required dependencies by following these steps:
- 2.1. Open the **Terminal**
- 2.2. Execute the following command:
```shell
git clone https://github.com/luanrabelo/SPLACE.git
cd SPLACE  
conda env create -f environment.yml
```
&nbsp;
- 2.3. Activate the conda environment:
```shell
conda activate splace_env
```
> [!NOTE]
> This will **clone the repository**, then you should navigate to the cloned directory to create the conda environment using the provided `environment.yml` file.
> After installation, activate the conda environment with: `conda activate splace_env`
&nbsp;  
### Usage
##### [:rocket: Go to Contents Overview](#contents-overview)
#### Parameter Overview
- `--input` : Path to the directory containing GenBank or Fasta files.
- `--output` : Path to the output directory where results will be saved. All marker files will be stored in zipped format, containing aligned and trimmed sequences, and a concatenated matrix, in NEXUS format for phylogenetic analysis.
- `-t` : Number of threads to use for parallel processing.
- `--gb_files` : Flag to indicate that the input files are in GenBank format.
- `--fasta_files` : Flag to indicate that the input files are in Fasta format.
- `--mafft` : Flag to enable multiple sequence alignment using MAFFT.
- `--trimal` : Flag to enable automated alignment trimming using TrimAl.
&nbsp;
#### Example Command
After installing **SPLACE** and activating the conda environment, you can run the tool using the command line interface. Here’s a basic example of how to use **SPLACE**:
- Run the following command in your terminal for processing **GenBank** files with **MAFFT** and **TrimAl**:
```shell
splace.py --input /path/to/your/gb_files --output /path/to/your/your_output_directory -t 8 --gb_files --mafft --trimal
```
- For processing **Fasta** files, use the following command:
```shell
splace.py --input /path/to/your/fasta_files --output /path/to/your/your_output_directory -t 8 --fasta_files --mafft --trimal
```
> [!NOTE]
> Make sure to replace `/path/to/your/gb_files`, `/path/to/your/fasta_files`, and `/path/to/your/your_output_directory` with the actual paths on your system.
> The `-t` parameter specifies the number of threads to use for parallel processing. Adjust this based on your system's capabilities.
> The `--mafft` and `--trimal` flags enable multiple sequence alignment and trimming, respectively. You can omit these flags if you do not wish to perform these steps.
> Do not use both `--gb_files` and `--fasta_files` flags together; choose one based on your input file format.
> For more detailed usage instructions and additional options, refer to the help command:
```shell
splace.py --help
```
&nbsp;
> [!CAUTION]
> For **Fasta files**, ensure that the sequence headers are formatted correctly to include gene names for proper processing by **SPLACE**.
> Example header format:
> ```
> > lcl|PX070005.1_cds_XZP64796.1_3 [gene=COX1] [protein=cytochrome c oxidase subunit I] [protein_id=XZP64796.1] [location=5509..7059] [gbkey=CDS]
> ATGGC...
> ```
> In this example, the gene name is specified as `COX1` ([gene=COX1]) and after standardization by **SynGenes**, it will be recognized as `COI` and written accordingly in the output files. This example header is a standard format for GenBank-derived Fasta files. 
> **Recommended practice is to use Fasta files generated from GenBank files to ensure compatibility.**
> Alternatively, for custom Fasta files, ensure the headers follow a similar structure to include the gene name, such as:
> ```
> > atp6_ITV1046I2 atp6 ATP synthase F0 subunit 6 7964:8638 forward
> ATGGC...
> ```
> In this case, the header splits into parts, using spaces as delimiters, where the second part (`atp6`) indicates the gene name. After standardization, it will be recognized as `ATP6` in the output files.
> If your Fasta headers, after splitting by spaces, do not have the gene name in the second position, SPALCE use the first part of the header as the gene name. This may lead to inconsistencies if the gene names are not standardized.
&nbsp;