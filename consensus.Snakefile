__author__ = "Taavi Päll"
__copyright__ = "Copyright 2020, Avilab"
__email__ = "taavi.pall@ut.ee"
__license__ = "MIT"


from snakemake.utils import validate, makedirs


# Load configuration file
configfile: "config.yaml"
validate(config, "schemas/config.schema.yaml")


REF_GENOME = config["refgenome"]
WRAPPER_PREFIX = "https://raw.githubusercontent.com/avilab/virome-wrappers/"


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
        expand("output/merged-{sample}/consensus.fa", sample = SAMPLES.keys())


rule vcfcombine:
    input:
        lambda wildcards: expand("output/{run}/lofreq.vcf", run = SAMPLES[wildcards.sample])
    output:
        "output/merged-{sample}/combined.vcf"
    resources:
        runtime = 120,
        mem_mb = 4000
    wrapper:
        WRAPPER_PREFIX + "master/vcflib/vcfcombine"


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
        WRAPPER_PREFIX + "master/vcflib/vcffilter"


rule referencemaker:
    input:
        vcf = "output/merged-{sample}/filtered.vcf",
        ref = REF_GENOME
    output:
        idx = temp("output/merged-{sample}/filtered.vcf.idx"),
        fasta = "output/merged-{sample}/consensus.fa",
        dic = "output/merged-{sample}/consensus.dict",
        fai = "output/merged-{sample}/consensus.fa.fai"
    params:
        refmaker = "--lenient"
    resources:
        runtime = 120,
        mem_mb = 4000    
    wrapper:
        WRAPPER_PREFIX + "master/gatk/fastaalternatereferencemaker"

rule rename:
    input:
        rules.referencemaker.output.fasta
    output:
        "output/merged-{sample}/consensus_rename.fa"
    params:
        sample = lambda wildcards: wildcards.sample,
        stub = "SARS-CoV-2/human/Estonia/"
    resources:
        runtime = 120,
        mem_mb = 2000    
    wrapper:
        WRAPPER_PREFIX + "master/sequences/rename_fasta"

rule merge_renamed:
    input:
        expand("output/merged-{sample}/consensus_rename.fa", sample = sample = SAMPLES.keys())
    output:
        "output/consensus.fa"
    resources:
        runtime = 120,
        mem_mb = 2000   
    shell:
        "cat {input} > {output}"
