__author__ = "Taavi Päll"
__copyright__ = "Copyright 2020, Avilab"
__email__ = "taavi.pall@ut.ee"
__license__ = "MIT"

# Load libraries
import os
import json
import glob
import pandas as pd
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
from snakemake.utils import validate, makedirs


# Load configuration file with sample and path info
configfile: "config.yaml"
validate(config, "schemas/config.schema.yaml")


# Load runs and groups
SAMPLES = pd.read_csv(config["samples"], sep="\s+")
validate(SAMPLES, "schemas/samples.schema.yaml")
SAMPLES = SAMPLES.set_index(
    ["run"], drop=False
)
RUN = SAMPLES.index.tolist()


# Path to reference genomes
REF_GENOME = config["refgenome"]
TAXON_DB = os.getenv("TAXON_DB")


# Wrappers
WRAPPER_PREFIX = "https://raw.githubusercontent.com/avilab/virome-wrappers/"


onsuccess:
    email = config["email"]
    shell("mail -s 'Forkflow finished successfully' {email} < {log}")


rule all:
    input: expand(["output/{run}/report.html", "output/{run}/genomecov.bg", "output/{run}/freebayes.vcf"], run = RUN)


def get_fastq(wildcards):
    paths = list(SAMPLES.loc[wildcards.run, ["fq1", "fq2"]])
    if config["remote"]:
        FTP = FTPRemoteProvider(username="anonymous", password=config["email"])
        return FTP.remote(paths, immediate_close=True)
    else:
        return paths


rule preprocess:
    input:
      sample = get_fastq
    output:
      adapters = temp("output/{run}/adapters.fa"),
      merged = temp("output/{run}/merged.fq"),
      unmerged = temp("output/{run}/unmerged.fq"),
      trimmed = temp("output/{run}/trimmed.fq"),
      sampled = temp("output/{run}/sample.fq")
    params:
      bbduk = "qtrim=r trimq=10 maq=10 minlen=100",
      seed = config["seed"]
    resources:
      runtime = 30,
      mem_mb = 8000
    threads: 4
    wrapper:
      WRAPPER_PREFIX + "master/preprocess"


# Map reads to ref genome
rule refgenome:
    input:
      reads = [rules.preprocess.output.sampled]
    output:
      "output/{run}/refgenome.bam"
    params:
      db_prefix = REF_GENOME,
      extra = "-L 100,100 -k 15",
      sorting = "samtools"
    resources:
      runtime = 30,
      mem_mb = 16000
    threads: 4
    wrapper:
      "https://raw.githubusercontent.com/tpall/snakemake-wrappers/bug/snakemake_issue145/bio/bwa/mem"


rule genomecov:
    input:
        ibam = rules.refgenome.output
    output:
        "output/{run}/genomecov.bg"
    params:
        extra = "-bg"
    resources:
      runtime = 20,
      mem_mb = 16000
    wrapper: 
        "file:../wrappers/bedtools/genomecov"


# Host mapping stats.
rule bamstats:
    input:
      rules.refgenome.output
    output:
      "output/{run}/bamstats.txt"
    params:
      extra = "-F 4",
      region = ""
    resources:
      runtime = 20,
      mem_mb = 8000
    wrapper:
      "0.42.0/bio/samtools/stats"


rule bcftools:
    input:
      ref=REF_GENOME,
      samples=rules.refgenome.output
    output:
      "output/{run}/bcftools.vcf"
    resources:
      runtime = 20,
      mem_mb = 8000
    params:
      mpileup = "-Ou --min-MQ 60",
      call = "-Ou -mv",
      norm = "-Ou -d all",
      filter = """-s LowQual -e '%QUAL<20'"""
    wrapper:
      "file:../wrappers/bcftools"


rule freebayes:
    input:
      ref=REF_GENOME,
      samples=rules.refgenome.output
    output:
        "output/{run}/freebayes.vcf" 
    params:
      extra="--ploidy 1"
    threads: 1
    wrapper:
        "file:../wrappers/freebayes"


rule report:
    input:
      bamstats = "output/{run}/bamstats.txt",
      vcf = "output/{run}/bcftools.vcf"
    output:
      "output/{run}/report.html"
    params:
      author = config["author"],
      run = lambda wildcards: wildcards.run
    resources:
      runtime = 10,
      mem_mb = 4000
    wrapper:
      "file:../wrappers/report"

