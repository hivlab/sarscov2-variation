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
df = pd.read_csv("samples.tsv", sep="\s+", dtype=str).set_index(["sample","run"], drop=False)
validate(df, "schemas/samples.schema.yaml")
samples = df.groupby(level=0).apply(lambda df: df.xs(df.name)["run"].tolist()).to_dict()
SAMPLE = [sample for sample,run in df.index.tolist()]
RUN = [run for sample,run in df.index.tolist()]
PLATFORM = "ILLUMINA"


# Path to reference genomes
REF_GENOME = config["refgenome"]
HOST_GENOME = os.environ["REF_GENOME_HUMAN_MASKED"]
RRNA_DB = os.environ["SILVA"]


# Wrappers
# Wrappers repo: https://github.com/avilab/virome-wrappers
WRAPPER_PREFIX = "https://raw.githubusercontent.com/avilab/virome-wrappers"


# Report
report: "report/workflow.rst"


onsuccess:
    email = config["email"]
    shell("mail -s 'Forkflow finished successfully' {email} < {log}")


rule all:
    input: 
        "output/consensus_masked.fa",
        "output/snpsift.csv", 
        "output/multiqc.html",
        expand(["output/{sample}/basecov.txt"], sample = list(samples.keys()))


def get_fastq(wildcards):
    fq_cols = [col for col in df.columns if "fq" in col]
    fqs = df.loc[(wildcards.sample, wildcards.run), fq_cols].dropna()
    assert len(fq_cols) in [1, 2], "Enter one or two FASTQ file paths"
    if len(fq_cols) == 2:
        return {"in1": fqs[0], "in2": fqs[1]}
    else:
        return {"input": fqs[0]}


rule reformat:
    input:
        unpack(get_fastq)
    output:
        out = temp("output/{sample}/{run}/interleaved.fq")
    params:
        extra = lambda wildcards, resources: f"-Xmx{resources.mem_mb / 1000:.0f}g da" 
    resources:
        runtime = 120,
        mem_mb = 4000
    log: 
        "output/{sample}/{run}/log/reformat.log"
    wrapper:
        f"{WRAPPER_PREFIX}/master/bbtools/reformat"


rule clumpify:
    input:
        input = rules.reformat.output.out
    output:
        out = temp("output/{sample}/{run}/clumpify.fq")
    params:
        extra = lambda wildcards, resources: f"-Xmx{resources.mem_mb / 1000:.0f}g dedupe optical spany adjacent markduplicates optical qin=33 da" # suppress assertions
    resources:
        runtime = 120,
        mem_mb = 4000
    threads: 8
    log: 
        "output/{sample}/{run}/log/clumpify.log"
    wrapper:
        f"{WRAPPER_PREFIX}/master/bbtools/clumpify"


rule trim:
    input:
        input = rules.clumpify.output.out
    output:
        out = temp("output/{sample}/{run}/trimmed.fq")
    params:
        extra = lambda wildcards, resources: f"-Xmx{resources.mem_mb / 1000:.0f}g maq=10 qtrim=r trimq=10 ktrim=r k=23 mink=11 hdist=1 tbo tpe minlen=100 ref=adapters,primers/primers.fa ftm=5 ordered qin=33"
    resources:
        runtime = 120,
        mem_mb = 4000
    log: 
        "output/{sample}/{run}/log/trim.log"
    wrapper:
        f"{WRAPPER_PREFIX}/master/bbtools/bbduk"


rule filter:
    input:
        input = rules.trim.output.out
    output:
        out = "output/{sample}/{run}/filtered.fq"
    params:
        extra = lambda wildcards, resources: f"-Xmx{resources.mem_mb / 1000:.0f}g k=31 ref=artifacts,phix ordered cardinality"
    resources:
        runtime = 120,
        mem_mb = 4000
    log: 
        "output/{sample}/{run}/log/filter.log"
    wrapper:
        f"{WRAPPER_PREFIX}/master/bbtools/bbduk"


# Remove rRNA sequences
rule maprRNA:
    input:
        input = rules.filter.output.out,
        ref = RRNA_DB
    output:
        outu = "output/{sample}/{run}/unmaprRNA.fq",
        outm = "output/{sample}/{run}/maprRNA.fq",
        statsfile = "output/{sample}/{run}/maprrna.txt"
    params:
        extra = lambda wildcards, resources: f"maxlen=600 nodisk -Xmx{resources.mem_mb / 1000:.0f}g"
    resources:
        runtime = 120,
        mem_mb = 16000
    threads: 4
    wrapper:
        f"{WRAPPER_PREFIX}/master/bbtools/bbwrap"


# Remove host sequences
rule maphost:
    input:
        input = rules.maprRNA.output.outu,
        ref = HOST_GENOME
    output:
        outu = "output/{sample}/{run}/unmaphost.fq",
        outm = "output/{sample}/{run}/maphost.fq",
        statsfile = "output/{sample}/{run}/maphost.txt"
    params:
        extra = lambda wildcards, resources: f"maxlen=600 nodisk -Xmx{resources.mem_mb / 1000:.0f}g"
    resources:
        runtime = lambda wildcards, attempt: attempt * 200,
        mem_mb = 24000
    threads: 4
    wrapper:
        f"{WRAPPER_PREFIX}/master/bbtools/bbwrap"


# Map reads to ref genome
rule refgenome:
    input:
        input = rules.maphost.output.outu,
        ref = REF_GENOME
    output:
        out = "output/{sample}/{run}/refgenome.bam",
        statsfile = "output/{sample}/{run}/refgenome.txt",
        gchist = "output/{sample}/{run}/gchist.txt",
        aqhist = "output/{sample}/{run}/aqhist.txt",
        lhist = "output/{sample}/{run}/lhist.txt",
        mhist = "output/{sample}/{run}/mhist.txt",
        bhist = "output/{sample}/{run}/bhist.txt",
    params:
        extra = lambda wildcards, resources: f"maxindel=200 strictmaxindel minid=0.9 maxlen=600 nodisk -Xmx{resources.mem_mb / 1000:.0f}g RGLB=lib1 RGPL={PLATFORM} RGID={wildcards.run} RGSM={wildcards.sample}"
    resources:
        runtime = 120,
        mem_mb = 16000
    threads: 4
    wrapper:
        f"{WRAPPER_PREFIX}/master/bbtools/bbwrap"


rule sort_and_index:
    input:
        rules.refgenome.output.out
    output:
        sorted = "output/{sample}/{run}/refgenome_sorted.bam",
        index = "output/{sample}/{run}/refgenome_sorted.bam.bai" 
    params:
        lambda wildcards, resources: f"-m {resources.mem_mb}M"
    threads:
        4
    resources:
        mem_mb = 16000,
        runtime = lambda wildcards, attempt: attempt * 240
    wrapper:
        f"{WRAPPER_PREFIX}/master/samtools/sort_and_index"


rule samtools_merge:
    input:
        lambda wildcards: expand("output/{{sample}}/{run}/refgenome_sorted.bam", run = samples[wildcards.sample])
    output:
        "output/{sample}/merged.bam"
    params:
        ""
    threads:  8  
    wrapper:
        "0.62.0/bio/samtools/merge"


rule pileup:
    input:
        input = rules.samtools_merge.output[0],
        ref = REF_GENOME
    output:
        out = "output/{sample}/covstats.txt",
        basecov = "output/{sample}/basecov.txt"
    params:
        extra = lambda wildcards, resources: f"-Xmx{resources.mem_mb}m concise"
    resources:
        runtime = lambda wildcards, attempt: attempt * 120,
        mem_mb = lambda wildcards, attempt: attempt * 8000
    wrapper: 
        f"{WRAPPER_PREFIX}/master/bbtools/pileup"


# Variant calling
# "Removes any sites with estimated probability of not being polymorphic 
# less than phred 20 (aka 0.01), or probability of polymorphism > 0.99"
# from FreeBayes user manual.
rule lofreq:
    input:
        ref = REF_GENOME,
        bam = rules.samtools_merge.output[0]
    output:
        "output/{sample}/lofreq.vcf" 
    params:
        extra="--call-indels --min-cov 50 --max-depth 1000000 --min-bq 30 --min-alt-bq 30 --def-alt-bq 0 --min-mq 20 --max-mq 255 --min-jq 0 --min-alt-jq 0 --def-alt-jq 0 --sig 0.01 --bonf dynamic --no-default-filter"
    resources:
        runtime = 120,
        mem_mb = 4000
    threads: 1
    wrapper:
        f"{WRAPPER_PREFIX}/master/lofreq/call"


rule vcffilter:
    input:
        "output/{sample}/lofreq.vcf"
    output:
        "output/{sample}/filtered.vcf"
    params:
        extra = "-f 'QUAL > 30 & AF > 0.5'"
    resources:
        runtime = 120,
        mem_mb = 4000
    wrapper:
        f"{WRAPPER_PREFIX}/master/vcflib/vcffilter"


rule genome_consensus:
    input:
        ref = REF_GENOME,
        bam = "output/{sample}/merged.bam",
        vcf = "output/{sample}/filtered.vcf"
    output:
        vcfgz = "output/{sample}/filtered.vcf.gz",
        consensus = "output/{sample}/consensus_badname.fa",
        consensus_masked = "output/{sample}/consensus_masked_badname.fa",
        bed = "output/{sample}/merged.bed"
    log:
        "output/{sample}/log/genome_consensus.log"
    params:
        mask = 20,
    wrapper:
        f"{WRAPPER_PREFIX}/master/genome-consensus"


rule rename:
    input:
        rules.genome_consensus.output.consensus_masked
    output:
        "output/{sample}/consensus_masked.fa"
    params:
        sample = lambda wildcards: wildcards.sample,
        stub = "SARS-CoV-2/human/Estonia/{}/2020"
    resources:
        runtime = 120,
        mem_mb = 2000    
    wrapper:
        f"{WRAPPER_PREFIX}/master/sequences/rename_fasta"


rule merge_renamed:
    input:
        expand("output/{sample}/consensus_masked.fa", sample = samples.keys())
    output:
        "output/consensus_masked.fa"
    resources:
        runtime = 120,
        mem_mb = 2000   
    shell:
        "cat {input} > {output}"


rule snpeff:
    input:
        "output/{sample}/lofreq.vcf"
    output:
        calls = "output/{sample}/snpeff.vcf",   # annotated calls (vcf, bcf, or vcf.gz)
        stats = "output/{sample}/snpeff.html",  # summary statistics (in HTML), optional
        csvstats = "output/{sample}/snpeff.csv", # summary statistics in CSV, optional
        genes = "output/{sample}/snpeff.genes.txt"
    log:
        "output/{sample}/log/snpeff_lofreq.log"
    params:
        data_dir = "data",
        reference = "NC045512", # reference name (from `snpeff databases`)
        extra = lambda wildcards, resources: f"-c refseq/snpEffect.config -Xmx{resources.mem_mb / 1000:.0f}g"          # optional parameters (e.g., max memory 4g)
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
        "output/{sample}/snpsift.txt"
    params:
        extra = "-s ',' -e '.'",
        fieldnames = "CHROM POS REF ALT DP AF SB DP4 EFF[*].IMPACT EFF[*].FUNCLASS EFF[*].EFFECT EFF[*].GENE EFF[*].CODON"
    wrapper:
        f"{WRAPPER_PREFIX}/master/snpsift"


rule merge_tables:
    input:
        expand("output/{sample}/snpsift.txt", sample = samples.keys())
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


# QC
fastq_screen_config = {
    "database": {
        "human": HOST_GENOME,
        "SILVA_138_SSU_132_LSU": RRNA_DB
    }
}


rule fastq_screen:
    input:
        rules.filter.output.out
    output:
        txt = "output/{sample}/{run}/fastq_screen.txt",
        png = "output/{sample}/{run}/fastq_screen.png"
    params:
        fastq_screen_config = fastq_screen_config,
        subset = 100000
    resources:
        runtime = 120,
        mem_mb = 8000    
    threads: 4
    wrapper:
        f"{WRAPPER_PREFIX}/master/fastq_screen"


rule fastqc:
    input:
        unpack(get_fastq)
    output:
        html = "output/{sample}/{run}/fastqc.html",
        zip = "output/{sample}/{run}/fastqc.zip"
    resources:
        runtime = 120,
        mem_mb = 4000    
    wrapper:
        "0.27.1/bio/fastqc"


# Host mapping stats
rule bamstats:
    input:
        rules.samtools_merge.output[0]
    output:
        "output/{sample}/bamstats.txt"
    resources:
        runtime = 120,
        mem_mb = 8000
    wrapper:
        "0.42.0/bio/samtools/stats"


rule multiqc:
    input:
        expand(["output/{sample}/{run}/fastq_screen.txt", "output/{sample}/{run}/fastqc.zip",], zip, sample = SAMPLE, run = RUN),
        expand(["output/{sample}/snpeff.csv", "output/{sample}/bamstats.txt"], sample = samples.keys())
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
      f"{WRAPPER_PREFIX}/master/multiqc"
