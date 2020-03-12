![status](https://img.shields.io/badge/status-under%20development-yellow)

# sarscov2

Snakemake workflow to align PE sequencing reads to NCBI reference sequence NC_045512.2.
Returns basic alignment stats and variants relative to refseq in html format.

## Installing

- Download and install miniconda3: <https://docs.conda.io/en/latest/miniconda.html>.
- Create conda environment and install snakemake.
```bash
conda create -n snakemake-env python=3.7
conda activate snakemake-env
conda install -c bioconda -c conda-forge snakemake
```

- Create a new working directory e.g. `covid-seq` and (Fork +) clone this repository to working directory.
```bash
mkdir covid-seq
cd covid-seq
git clone https://github.com/avilab/sarscov2.git .
```

- Edit `samples.tsv` with full paths to sequencing reads and run names and edit also `config.yaml`.

## Running
Test run:
```bash
snakemake --use-conda -n
```

Analyse sequences:
```bash
snakemake --use-conda -j
```

For all possible snakemake command line options please refer to snakemake tutorial <https://snakemake.readthedocs.io/en/stable/executing/cli.html>.


This workflow can be run on a contemporary PC/laptop (e.g. i5/16G) with sufficient HD space to accomodate sequening runs.



