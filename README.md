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

- Create a new directory e.g. `covid-seq` and (Fork +) clone this repository.
```bash
mkdir covid-seq
cd covid-seq
git clone avilab/sarscov2 
```
