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
SAMPLES = pd.read_csv(config["samples"], sep="\s+", dtype=str).set_index(["run"], drop=False)
validate(SAMPLES, "schemas/samples.schema.yaml")
RUN = SAMPLES.index.tolist()


# Path to reference genomes
REF_GENOME = config["refgenome"]
HOST_GENOME = os.environ["REF_GENOME_HUMAN_MASKED"]
RRNA_DB = os.environ["SILVA"]


# Wrappers
WRAPPER_PREFIX = "https://raw.githubusercontent.com/avilab/virome-wrappers/"


# Report
report: "report/workflow.rst"


onsuccess:
    email = config["email"]
    shell("mail -s 'Forkflow finished successfully' {email} < {log}")


rule all:
    input: expand(["output/{run}/filtered.fq", "output/{run}/unmaphost.fq", "output/{run}/fastq_screen.txt", "output/{run}/fastqc.zip"], run = RUN)


def get_fastq(wildcards):
    fq_cols = [col for col in SAMPLES.columns if "fq" in col]
    fqs = SAMPLES.loc[wildcards.run, fq_cols].dropna()
    assert len(fq_cols) in [1, 2], "Enter one or two FASTQ file paths"
    if len(fq_cols) == 2:
        return {"in1": fqs[0], "in2": fqs[1]}
    else:
        return {"input": fqs[0]}


rule clumpify:
    input:
        unpack(get_fastq)
    output:
        out = temp("output/{run}/clumpify.fq")
    params:
        extra = "dedupe optical qin=33 -da" # suppress assertions
    resources:
        runtime = 20,
        mem_mb = 4000
    log: 
        "output/{run}/log/clumpify.log"
    wrapper:
        WRAPPER_PREFIX + "master/bbmap/clumpify"


rule trim:
    input:
        input = rules.clumpify.output.out
    output:
        out = temp("output/{run}/trimmed.fq")
    params:
        extra = "ktrim=r k=23 mink=11 hdist=1 tbo tpe minlen=70 ref=adapters ftm=5 ordered qin=33"
    resources:
        runtime = 20,
        mem_mb = 4000
    log: 
        "output/{run}/log/trim.log"
    wrapper:
        WRAPPER_PREFIX + "master/bbduk"


rule filter:
    input:
        input = rules.trim.output.out
    output:
        out = "output/{run}/filtered.fq"
    params:
        extra = "k=31 ref=artifacts,phix ordered cardinality"
    resources:
        runtime = 20,
        mem_mb = 4000
    log: 
        "output/{run}/log/filter.log"
    wrapper:
        WRAPPER_PREFIX + "master/bbduk"


# Remove rRNA sequences
rule maprRNA:
    input:
        input = rules.filter.output.out,
        ref = RRNA_DB
    output:
        outu = "output/{run}/unmaprRNA.fq",
        outm = "output/{run}/maprRNA.fq",
        statsfile = "output/{run}/maprrna.txt"
    params:
        extra = "maxlen=600 nodisk -Xmx16000m"
    resources:
        runtime = 30,
        mem_mb = 16000
    threads: 4
    wrapper:
        WRAPPER_PREFIX + "master/bbmap/bbwrap"


# Remove host sequences
rule maphost:
    input:
        input = rules.maprRNA.output.outu,
        ref = HOST_GENOME
    output:
        outu = "output/{run}/unmaphost.fq",
        outm = "output/{run}/maphost.fq",
        statsfile = "output/{run}/maphost.txt"
    params:
        extra = "maxlen=600 nodisk -Xmx16000m"
    resources:
        runtime = 30,
        mem_mb = 16000
    threads: 4
    wrapper:
        WRAPPER_PREFIX + "master/bbmap/bbwrap"


# QC
fastq_screen_config = {
    "database": {
        "human": HOST_GENOME,
        "SILVA_138_SSU_132_LSU": RRNA_DB
    }
}
rule fastq_screen:
    input:
        rules.trim.output.out
    output:
        txt = "output/{run}/fastq_screen.txt",
        png = "output/{run}/fastq_screen.png"
    params:
        fastq_screen_config = fastq_screen_config,
        subset = 100000
    resources:
        runtime = 30,
        mem_mb = 8000    
    threads: 4
    wrapper:
        WRAPPER_PREFIX + "master/fastq_screen"


rule fastqc:
    input:
        unpack(get_fastq)
    output:
        html = "output/{run}/fastqc.html",
        zip = "output/{run}/fastqc.zip"
    resources:
        runtime = 20,
        mem_mb = 4000    
    wrapper:
        "0.27.1/bio/fastqc"

