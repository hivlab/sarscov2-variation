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
    input: expand(["output/{run}/report.html", "output/{run}/final.contigs.fa", "output/{run}/coverage.txt", "output/{run}/basecov.txt", "output/{run}/genomecov.bg"], run = RUN)



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
    threads: 4
    wrapper:
      WRAPPER_PREFIX + "master/preprocess"


# Map reads to ref genome
rule bwa_mem_ref:
    input:
      reads = [rules.preprocess.output.sampled]
    output:
      "output/{run}/refgenome.bam"
    params:
      db_prefix = REF_GENOME,
      extra = "-L 100,100 -k 15",
      sorting = "samtools"
    threads: 4
    wrapper:
      "https://raw.githubusercontent.com/tpall/snakemake-wrappers/bug/snakemake_issue145/bio/bwa/mem"


rule genomecov:
    input:
        ibam = rules.bwa_mem_ref.output
    output:
        "output/{run}/genomecov.bg"
    params:
        extra = "-bg"
    wrapper: 
        "file:../wrappers/bedtools/genomecov"


rule ref_mapped:
    input:
      rules.bwa_mem_ref.output
    output:
      "output/{run}/mapped.bam"
    params:
      "-b -F 4"
    threads: 4
    wrapper:
      "0.49.0/bio/samtools/view"


# Host mapping stats.
rule ref_bam_stats:
    input:
      rules.bwa_mem_ref.output
    output:
      "output/{run}/bamstats.txt"
    params:
      extra = "-F 4",
      region = ""
    wrapper:
      "0.42.0/bio/samtools/stats"


rule bcftools_call:
    input:
      ref=REF_GENOME,
      samples=rules.bwa_mem_ref.output
    output:
      "output/{run}/var_filt.vcf"
    params:
      mpileup="-Ou",
      call="-Ou -mv",
      filter="-s LowQual -e '%QUAL<20 || DP>100'"
    wrapper:
      "file:../wrappers/bcftools"


rule samtools_bam2fq:
    input:
      rules.ref_mapped.output
    output:
      "output/{run}/mapped.fq"
    threads: 4
    wrapper:
        "0.49.0/bio/samtools/bam2fq/interleaved"


rule assemble:
    input: 
      se = expand("output/{run}/mapped.fq", run = RUN)
    output: 
      contigs = "output/{run}/final.contigs.fa"
    params:
      extra = ""
    threads: 4
    log: "logs/{run}_assemble.log"
    shadow: 
      "minimal"
    wrapper:
      WRAPPER_PREFIX + "release/metformin-pill/assembly/megahit"


# Calculate assembly coverage stats
# nodisk keeps index in memory, otherwise index will be written once to project root (ref/1) from first run to be processed 
# and reused for other unrelated runs.
# Key "input" will be parsed to "in", "input1" to "in1" etc.
rule coverage:
    input:
      ref = rules.assemble.output.contigs, 
      input = expand("output/{run}/mapped.fq", run = RUN) 
    output:
      out = temp("output/{run}/final.contigs_aln.sam"),
      covstats = "output/{run}/coverage.txt",
      basecov = "output/{run}/basecov.txt"
    params: 
      extra = "nodisk"
    wrapper:
      WRAPPER_PREFIX + "master/bbmap/bbwrap"


rule report:
    input:
      bamstats = "output/{run}/bamstats.txt",
      vcf = "output/{run}/var_filt.vcf"
    output:
      "output/{run}/report.html"
    params:
      author = config["author"],
      run = lambda wildcards: wildcards.run
    wrapper:
      "file:../wrappers/report"

