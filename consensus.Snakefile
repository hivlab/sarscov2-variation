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
CALLER = ["freebayes", "lofreq"]


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
        expand("output/{sample}/{caller}_consensus.fa", sample = SAMPLES.keys(), caller = CALLER)


rule vcfcombine:
    input:
        lambda wildcards: expand("output/{run}/{{caller}}.vcf", run = SAMPLES[wildcards.sample])
    output:
        "output/{sample}/{caller}_combined.vcf"
    resources:
        runtime = 120,
        mem_mb = 4000
    wrapper:
        WRAPPER_PREFIX + "master/vcflib/vcfcombine"


rule vcffilter:
    input:
        "output/{sample}/{caller}_combined.vcf"
    output:
        "output/{sample}/{caller}_filtered.vcf"
    params:
        extra = "-f 'QUAL > 30 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1'"
    resources:
        runtime = 120,
        mem_mb = 4000
    wrapper:
        WRAPPER_PREFIX + "master/vcflib/vcffilter"


rule referencemaker:
    input:
        vcf = "output/{sample}/{caller}_combined.vcf",
        ref = REF_GENOME
    output:
        idx = temp("output/{sample}/{caller}_combined.vcf.idx"),
        fasta = "output/{sample}/{caller}_consensus.fa",
        dic = "output/{sample}/{caller}_consensus.dict",
        fai = "output/{sample}/{caller}_consensus.fa.fai"
    params:
        refmaker = "--lenient"
    resources:
        runtime = 120,
        mem_mb = 4000    
    wrapper:
        WRAPPER_PREFIX + "master/gatk/fastaalternatereferencemaker"

