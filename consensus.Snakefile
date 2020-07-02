__author__ = "Taavi Päll"
__copyright__ = "Copyright 2020, Avilab"
__email__ = "taavi.pall@ut.ee"
__license__ = "MIT"


from snakemake.utils import validate, makedirs


# Load configuration file
configfile: "config.yaml"
validate(config, "schemas/config.schema.yaml")


REF_GENOME = config["refgenome"]
WRAPPER_PREFIX = "https://raw.githubusercontent.com/avilab/virome-wrappers"


SAMPLES = {
    "LB202": ["LB202", "LB202-I"],
    "LB614": ["LB614", "LB614-I"],
    "LB624": ["LB624", "LB624-I"],
    "LB228": ["LB228", "LB228-I"],
    "LB637": ["LB637", "LB637-I"],
    "LB236": ["LB236", "LB236-I"],
}


rule all:
    input:
        expand("output/merged-{sample}/consensus.fa", sample = SAMPLES.keys()), "output/consensus.fa"


rule vcfcombine:
    input:
        lambda wildcards: expand("output/{run}/lofreq.vcf", run = SAMPLES[wildcards.sample])
    output:
        "output/merged-{sample}/combined.vcf"
    resources:
        runtime = 120,
        mem_mb = 4000
    wrapper:
        f"{WRAPPER_PREFIX}/master/vcflib/vcfcombine"


rule vcffilter:
    input:
        "output/merged-{sample}/combined.vcf"
    output:
        "output/merged-{sample}/filtered.vcf"
    params:
        extra = "-f 'QUAL > 30 & AF > 0.5'"
    resources:
        runtime = 120,
        mem_mb = 4000
    wrapper:
        f"{WRAPPER_PREFIX}/master/vcflib/vcffilter"


rule samtools_merge:
    input:
        lambda wildcards: expand("output/{run}/refgenome_sorted.bam", run = SAMPLES[wildcards.sample])
    output:
        "output/merged-{sample}/merged.bam"
    params:
        ""
    threads:  8  
    wrapper:
        "0.62.0/bio/samtools/merge"


rule genome_consensus:
    input:
        ref = REF_GENOME,
        bam = "output/merged-{sample}/merged.bam",
        vcf = "output/merged-{sample}/filtered.vcf"
    output:
        vcfgz = "output/merged-{sample}/filtered.vcf.gz",
        consensus = "output/merged-{sample}/consensus_badname.fa",
        consensus_masked = "output/merged-{sample}/consensus_masked_badname.fa",
        bed = "output/merged-{sample}/merged.bed"
    log:
        "output/merged-{sample}/log/genome_consensus.log"
    params:
        mask = 20,
    wrapper:
        f"{WRAPPER_PREFIX}/master/bbtools/bbduk"


rule rename:
    input:
        rules.genome_consensus.output.consensus_masked
    output:
        "output/merged-{sample}/consensus_masked.fa"
    params:
        sample = lambda wildcards: wildcards.sample,
        stub = "SARS-CoV-2/human/Estonia/{}/2020"
    resources:
        runtime = 120,
        mem_mb = 2000    
    wrapper:
        "file:wrappers/sequences/rename_fasta"

rule merge_renamed:
    input:
        expand("output/merged-{sample}/consensus_masked.fa", sample = SAMPLES.keys())
    output:
        "output/consensus_masked.fa"
    resources:
        runtime = 120,
        mem_mb = 2000   
    shell:
        "cat {input} > {output}"
