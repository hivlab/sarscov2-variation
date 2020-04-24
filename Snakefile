__author__ = "Taavi Päll"
__copyright__ = "Copyright 2020, Avilab"
__email__ = "taavi.pall@ut.ee"
__license__ = "MIT"

# Load libraries
import os
import pandas as pd
from snakemake.utils import validate, makedirs


# Load configuration file with sample and path info
configfile: "config.yaml"
validate(config, "schemas/config.schema.yaml")


# Load runs and groups
SAMPLES = pd.read_csv(config["samples"], sep="\s+", dtype=str).set_index(["run"], drop=False)
validate(SAMPLES, "schemas/samples.schema.yaml")
RUN = SAMPLES.index.tolist()
PLATFORM = SAMPLES.platform


# Path to reference genomes
REF_GENOME = config["refgenome"]
HOST_GENOME = os.environ["REF_GENOME_HUMAN_MASKED"]
RRNA_DB = os.environ["SILVA"]
# cpn60 NR database file cpndb_nr_nut_seq.txt was downloaded from http://www.cpndb.ca/downloads.php
# cpn60 database was indexed using bwa index
CPNDB = os.environ["CPNDB"]


# Wrappers
# Wrappers repo: https://github.com/avilab/virome-wrappers
WRAPPER_PREFIX = "https://raw.githubusercontent.com/avilab/virome-wrappers/"


# Report
report: "report/workflow.rst"


onsuccess:
    email = config["email"]
    shell("mail -s 'Forkflow finished successfully' {email} < {log}")


rule all:
    input: expand(["output/snpsift_lofreq.csv", "output/snpsift.csv", "output/{run}/report.html", "output/multiqc.html", "output/{run}/freebayes.vcf", "output/{run}/filtered.fq", "output/{run}/unmaphost.fq", "output/{run}/fastq_screen.txt", "output/{run}/fastqc.zip"], run = RUN)

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
        extra = "dedupe optical reorder qin=33 -da" # suppress assertions
    resources:
        runtime = 120,
        mem_mb = 4000
    log: 
        "output/{run}/log/clumpify.log"
    wrapper:
        WRAPPER_PREFIX + "master/bbtools/clumpify"


rule trim:
    input:
        input = rules.clumpify.output.out
    output:
        out = temp("output/{run}/trimmed.fq")
    params:
        extra = "ktrim=r k=23 mink=11 hdist=1 tbo tpe minlen=70 ref=adapters,primers/primers.fa ftm=5 ordered qin=33"
    resources:
        runtime = 120,
        mem_mb = 4000
    log: 
        "output/{run}/log/trim.log"
    wrapper:
        WRAPPER_PREFIX + "master/bbtools/bbduk"


rule filter:
    input:
        input = rules.trim.output.out
    output:
        out = "output/{run}/filtered.fq"
    params:
        extra = "k=31 ref=artifacts,phix ordered cardinality"
    resources:
        runtime = 120,
        mem_mb = 4000
    log: 
        "output/{run}/log/filter.log"
    wrapper:
        WRAPPER_PREFIX + "master/bbtools/bbduk"


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
        extra = "maxlen=600 nodisk -Xmx16g"
    resources:
        runtime = 120,
        mem_mb = 16000
    threads: 4
    wrapper:
        WRAPPER_PREFIX + "master/bbtools/bbwrap"


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
        extra = "maxlen=600 nodisk -Xmx24g"
    resources:
        runtime = lambda wildcards, attempt: attempt * 200,
        mem_mb = 24000
    threads: 4
    wrapper:
        WRAPPER_PREFIX + "master/bbtools/bbwrap"


# Map reads to ref genome
rule refgenome:
    input:
        input = rules.maphost.output.outu,
        ref = REF_GENOME
    output:
        out = "output/{run}/refgenome.sam",
        statsfile = "output/{run}/refgenome.txt",
        gchist = "output/{run}/gchist.txt",
        aqhist = "output/{run}/aqhist.txt",
        lhist = "output/{run}/lhist.txt",
        mhist = "output/{run}/mhist.txt",
        bhist = "output/{run}/bhist.txt",
    params:
        extra = "maxlen=600 nodisk -Xmx16g"
    resources:
        runtime = 120,
        mem_mb = 16000
    threads: 4
    wrapper:
        WRAPPER_PREFIX + "master/bbtools/bbwrap"


rule replace_rg:
    input:
        rules.refgenome.output.out
    output:
        "output/{run}/refgenome_fixed.bam"
    params:
        lambda wildcards: "RGLB=lib1 RGPL={} RGPU={} RGSM={}".format(PLATFORM[wildcards.run], wildcards.run, wildcards.run)
    resources:
        runtime = 120,
        mem_mb = 4000
    wrapper:
        WRAPPER_PREFIX + "master/picard/addorreplacereadgroups"


rule samtools_sort:
    input:
        rules.replace_rg.output[0]
    output:
        "output/{run}/refgenome_sorted.bam"
    params:
        ""
    resources:
        runtime = 120,
        mem_mb = 4000
    threads: 4 # Samtools takes additional threads through its option -@
    wrapper:
        "0.50.4/bio/samtools/sort"

rule genomecov:
    input:
        ibam = rules.replace_rg.output[0]
    output:
        "output/{run}/genomecov.bg"
    params:
        extra = "-bg"
    resources:
        runtime = 120,
        mem_mb = lambda wildcards, attempt: attempt * 8000
    wrapper: 
        WRAPPER_PREFIX + "master/bedtools/genomecov"


# Variant calling
# "Removes any sites with estimated probability of not being polymorphic 
# less than phred 20 (aka 0.01), or probability of polymorphism > 0.99"
# from FreeBayes user manual.
rule freebayes:
    input:
        ref = REF_GENOME,
        samples = rules.samtools_sort.output[0]
    output:
        "output/{run}/freebayes.vcf" 
    params:
        extra="--pooled-continuous --ploidy 1",
        pipe = ""
    resources:
        runtime = 120,
        mem_mb = 4000
    threads: 1
    wrapper:
        WRAPPER_PREFIX + "master/freebayes"


rule lofreq:
    input:
        ref = REF_GENOME,
        bam = rules.samtools_sort.output[0]
    output:
        "output/{run}/lofreq.vcf" 
    params:
        extra="--call-indels --min-cov 50 --min-bq 30 --min-alt-bq 30 --no-ext-baq --min-mq 20 --max-mq 255"
    resources:
        runtime = 120,
        mem_mb = 4000
    threads: 1
    wrapper:
        WRAPPER_PREFIX + "master/lofreq"


rule snpeff:
    input:
        "output/{run}/freebayes.vcf"
    output:
        calls = "output/{run}/snpeff.vcf",   # annotated calls (vcf, bcf, or vcf.gz)
        stats = "output/{run}/snpeff.html",  # summary statistics (in HTML), optional
        csvstats = "output/{run}/snpeff.csv", # summary statistics in CSV, optional
        genes = "output/{run}/snpeff.genes.txt"
    log:
        "output/{run}/log/snpeff.log"
    params:
        data_dir = "data",
        reference = "NC045512", # reference name (from `snpeff databases`)
        extra = "-c refseq/snpEffect.config -Xmx4g"          # optional parameters (e.g., max memory 4g)
    resources:
        runtime = 120,
        mem_mb = 4000    
    wrapper:
        "0.50.4/bio/snpeff"


rule snpeff_lofreq:
    input:
        "output/{run}/lofreq.vcf"
    output:
        calls = "output/{run}/snpeff_lofreq.vcf",   # annotated calls (vcf, bcf, or vcf.gz)
        stats = "output/{run}/snpeff_lofreq.html",  # summary statistics (in HTML), optional
        csvstats = "output/{run}/snpeff_lofreq.csv", # summary statistics in CSV, optional
        genes = "output/{run}/snpeff_lofreq.genes.txt"
    log:
        "output/{run}/log/snpeff_lofreq.log"
    params:
        data_dir = "data",
        reference = "NC045512", # reference name (from `snpeff databases`)
        extra = "-c refseq/snpEffect.config -Xmx4g"          # optional parameters (e.g., max memory 4g)
    resources:
        runtime = 120,
        mem_mb = 4000    
    wrapper:
        "0.50.4/bio/snpeff"



# Parse snpeff output to tabular format
rule snpsift:
    input:
        rules.snpeff.output.calls
    output:
        "output/{run}/snpsift.txt"
    params:
        extra = "-s ',' -e '.'",
        fieldnames = "CHROM POS REF ALT QUAL AO SAF SAR RPR RPL FILTER AF SB DP ANN[*].IMPACT ANN[*].EFFECT ANN[*].GENE ANN[*].CODON"
    wrapper:
        WRAPPER_PREFIX + "master/snpsift"


rule snpsift_lofreq:
    input:
        rules.snpeff_lofreq.output.calls
    output:
        "output/{run}/snpsift_lofreq.txt"
    params:
        extra = "-s ',' -e '.'",
        fieldnames = "CHROM POS REF ALT DP AF SB DP4 EFF[*].IMPACT EFF[*].FUNCLASS EFF[*].EFFECT EFF[*].GENE EFF[*].CODON"
    wrapper:
        WRAPPER_PREFIX + "master/snpsift"


rule merge_tables:
    input:
        expand("output/{run}/snpsift.txt", run = RUN)
    output:
        "output/snpsift.csv"
    run:
        import pandas as pd
        files = {}
        for file in input:
            files.update({file.split("/")[1]: pd.read_csv(file, sep = "\t")})
        concatenated = pd.concat(files, names = ["Sample"])
        modified = concatenated.reset_index()
        modified.to_csv(output[0], index = False)


rule merge_tables_lofreq:
    input:
        expand("output/{run}/snpsift_lofreq.txt", run = RUN)
    output:
        "output/snpsift_lofreq.csv"
    run:
        import pandas as pd
        files = {}
        for file in input:
            files.update({file.split("/")[1]: pd.read_csv(file, sep = "\t")})
        concatenated = pd.concat(files, names = ["Sample"])
        modified = concatenated.reset_index()
        modified.to_csv(output[0], index = False)


# Parse report
rule report:
    input:
        statsfile = "output/{run}/refgenome.txt",
        gchist = "output/{run}/gchist.txt",
        aqhist = "output/{run}/aqhist.txt",
        lhist = "output/{run}/lhist.txt",
        mhist = "output/{run}/mhist.txt",
        bhist = "output/{run}/bhist.txt",
        genomecov = "output/{run}/genomecov.bg",
        vcf = "output/{run}/freebayes.vcf"
    output:
        "output/{run}/report.html"
    params:
        author = config["author"],
        run = lambda wildcards: wildcards.run
    log: 
        "output/{run}/log/report.log"
    resources:
        runtime = 120,
        mem_mb = 4000
    wrapper:
        "file:wrappers/report"

# QC
fastq_screen_config = {
    "database": {
        "human": HOST_GENOME,
        "SILVA_138_SSU_132_LSU": RRNA_DB,
        "cpn60": CPNDB
    }
}

rule fastq_screen:
    input:
        rules.filter.output.out
    output:
        txt = "output/{run}/fastq_screen.txt",
        png = "output/{run}/fastq_screen.png"
    params:
        fastq_screen_config = fastq_screen_config,
        subset = 100000
    resources:
        runtime = 120,
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
        runtime = 120,
        mem_mb = 4000    
    wrapper:
        "0.27.1/bio/fastqc"


# Host mapping stats
rule bamstats:
    input:
        rules.replace_rg.output
    output:
        "output/{run}/bamstats.txt"
    resources:
        runtime = 120,
        mem_mb = 8000
    wrapper:
        "0.42.0/bio/samtools/stats"


rule multiqc:
    input:
        expand(["output/{run}/fastq_screen.txt",
        "output/{run}/bamstats.txt",
        "output/{run}/fastqc.zip",
        "output/{run}/snpeff.csv"], run = RUN)
    output:
        report("output/multiqc.html", caption = "report/multiqc.rst", category = "Quality control")
    params:
        "-d -dd 1"
    log:
        "output/multiqc.log"
    resources:
        runtime = 120,
        mem_mb = 4000    
    wrapper:
      WRAPPER_PREFIX + "master/multiqc"
