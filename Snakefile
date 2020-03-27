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
    input: expand(["output/{run}/report.html", "output/{run}/genomecov.bg", "output/{run}/freebayes.vcf", "output/{run}/all.vcf", "output/{run}/sample.fq"], run = RUN)


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
      sampled = "output/{run}/sample.fq"
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
      input = [rules.preprocess.output.sampled],
      ref = REF_GENOME
    output:
      out = "output/{run}/refgenome.sam",
      statsfile = "output/{run}/statsfile.txt",
      bhist = "output/{run}/bhist.txt",
      qhist = "output/{run}/qhist.txt",
      aqhist = "output/{run}/aqhist.txt",
      lhist = "output/{run}/lhist.txt", 
      ihist = "output/{run}/ihist.txt", 
      ehist = "output/{run}/ehist.txt", 
      qahist = "output/{run}/qahist.txt", 
      indelhist = "output/{run}/indelhist.txt", 
      mhist = "output/{run}/mhist.txt", 
      gchist = "output/{run}/gchist.txt", 
      idhist = "output/{run}/idhist.txt"
    params:
      extra = "maxlen=600 nodisk"
    resources:
      runtime = 30,
      mem_mb = 16000
    threads: 4
    wrapper:
      WRAPPER_PREFIX + "master/bbmap/bbwrap"


rule samtools_sort:
    input:
      rules.refgenome.output.out
    output:
      "output/{run}/refgenome.bam"
    params:
      "-m 4G"
    threads: 4 # Samtools takes additional threads through its option -@
    wrapper:
        "0.50.4/bio/samtools/sort"


rule replace_rg:
    input:
      rules.samtools_sort.output
    output:
      "output/{run}/refgenome_fixed.bam"
    params:
      "RGLB=lib1 RGPL=ILLUMINA RGPU={run} RGSM={run}"
    wrapper:
      WRAPPER_PREFIX + "master/picard/addorreplacereadgroups"


rule dedup:
    input: 
      rules.replace_rg.output
    output:
      bam = "output/{run}/refgenome_dedup.bam",
      metrics = "output/{run}/dedup.metrics"
    params:
      extra = "REMOVE_DUPLICATES=TRUE ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT"
    wrapper:
      WRAPPER_PREFIX + "master/picard/markduplicates"


rule genomecov:
    input:
      ibam = rules.dedup.output.bam
    output:
      "output/{run}/genomecov.bg"
    params:
      extra = "-bg"
    resources:
      runtime = 20,
      mem_mb = 16000
    wrapper: 
      WRAPPER_PREFIX + "master/bedtools/genomecov"


# Host mapping stats.
rule bamstats:
    input:
      rules.dedup.output.bam
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
      samples=rules.dedup.output.bam
    output:
      "output/{run}/bcftools.vcf"
    resources:
      runtime = 20,
      mem_mb = 8000
    params:
      mpileup = "-Ou --min-MQ 60",
      call = "-Ou -mv --ploidy 1",
      norm = "-Ou -d all",
      filter = """-s LowQual -e '%QUAL<20'"""
    wrapper:
      WRAPPER_PREFIX + "master/bcftools"


rule freebayes:
    input:
      ref=REF_GENOME,
      samples=rules.dedup.output.bam
    output:
        "output/{run}/freebayes.vcf" 
    params:
      extra="--ploidy 1",
      pipe = """| bcftools filter -s LowQual -e '%QUAL<20' """
    threads: 1
    wrapper:
      WRAPPER_PREFIX + "master/freebayes"


rule bcftools_concat:
    input:
        calls=["output/{run}/bcftools.vcf", "output/{run}/freebayes.vcf"]
    output:
        "output/{run}/all.vcf"
    params:
        ""  # optional parameters for bcftools concat (except -o)
    wrapper:
        "0.50.4/bio/bcftools/concat"


rule report:
    input:
      bamstats = "output/{run}/bamstats.txt",
      vcf = "output/{run}/all.vcf"
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

