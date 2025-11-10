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
    - `syngenes`  # For gene nomenclature standardization
    - `biopython` # For biological sequence handling and parsing
    - `mafft`     # For multiple sequence alignment
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
> [!NOTE]
> This will **clone the repository**, then you should navigate to the cloned directory to create the conda environment using the provided `environment.yml` file.
> After installation, activate the conda environment with: `conda activate splace_env`
&nbsp;  