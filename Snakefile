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
HOST_GENOME = os.environ["REF_GENOME_HUMAN_MASKED"]
RRNA_DB = os.environ["SILVA_SSU"]


# Wrappers
WRAPPER_PREFIX = "https://raw.githubusercontent.com/avilab/virome-wrappers/"


onsuccess:
    email = config["email"]
    shell("mail -s 'Forkflow finished successfully' {email} < {log}")


rule all:
    input: expand(["output/{run}/report.html", "output/{run}/genomecov.bg", "output/{run}/freebayes.vcf", "output/{run}/all.vcf", "output/{run}/rrna_statsfile.txt", "output/{run}/host_statsfile.txt"], run = RUN)


def get_fastq(wildcards):
    paths = list(SAMPLES.loc[wildcards.run, ["fq1", "fq2"]])
    if config["remote"]:
        FTP = FTPRemoteProvider(username="anonymous", password=config["email"])
        return FTP.remote(paths, immediate_close=True)
    else:
        return paths


rule preprocess:
    input:
      in1 = lambda wildcards: get_fastq(wildcards)[0],
      in2 = lambda wildcards: get_fastq(wildcards)[1]
    output:
      out1 = temp("output/{run}/cleaned1.fq"),
      out2 = temp("output/{run}/cleaned2.fq"),
      gchist = "output/{run}/pre_gchist.txt",
      aqhist = "output/{run}/pre_aqhist.txt",
      lhist = "output/{run}/pre_lhist.txt",
      mhist = "output/{run}/pre_mhist.txt",
      bhist = "output/{run}/pre_bhist.txt"
    params:
      extra = "hdist=1 maq=20 tpe tbo -da"
    resources:
      runtime = 20,
      mem_mb = 4000
    log: "output/{run}/bbduk.log"
    wrapper:
      WRAPPER_PREFIX + "master/bbduk"


# Map reads to ref genome
rule refgenome:
    input:
      in1 = rules.preprocess.output.out1,
      in2 = rules.preprocess.output.out2,
      ref = REF_GENOME
    output:
      out = "output/{run}/refgenome.sam",
      statsfile = "output/{run}/statsfile.txt",
      bhist = "output/{run}/bhist.txt",
      aqhist = "output/{run}/aqhist.txt",
      lhist = "output/{run}/lhist.txt", 
      mhist = "output/{run}/mhist.txt", 
      gchist = "output/{run}/gchist.txt"
    params:
      extra = "maxlen=600 nodisk -Xmx8000m"
    resources:
      runtime = 30,
      mem_mb = 8000
    threads: 4
    wrapper:
      WRAPPER_PREFIX + "master/bbmap/bbwrap"


rule rrna:
    input:
      in1 = rules.preprocess.output.out1,
      in2 = rules.preprocess.output.out2,
      ref = RRNA_DB
    output:
      out = temp("output/{run}/norrna.sam"),
      statsfile = "output/{run}/rrna_statsfile.txt"
    params:
      extra = "maxlen=600 nodisk -Xmx16000m"
    resources:
      runtime = 30,
      mem_mb = 16000
    threads: 4
    wrapper:
      WRAPPER_PREFIX + "master/bbmap/bbwrap"


rule host:
    input:
      in1 = rules.preprocess.output.out1,
      in2 = rules.preprocess.output.out2,
      ref = HOST_GENOME
    output:
      out = temp("output/{run}/host.sam"),
      statsfile = "output/{run}/host_statsfile.txt"
    params:
      extra = "maxlen=600 nodisk -Xmx16000m"
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


# Host mapping stats
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


# Variant calling
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


# Parse report
rule report:
    input:
      statsfile = "output/{run}/statsfile.txt",
      gchist = "output/{run}/gchist.txt",
      aqhist = "output/{run}/aqhist.txt",
      lhist = "output/{run}/lhist.txt",
      mhist = "output/{run}/mhist.txt",
      bhist = "output/{run}/bhist.txt",
      genomecov = "output/{run}/genomecov.bg",
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

