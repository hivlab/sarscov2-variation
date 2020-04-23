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
    "LB775": ["LB775"],
    "LC039": ["LC039"],
}

rule all:
    input:
        expand("output/{sample}/consensus.fa", sample = SAMPLES.keys())


rule vcfcombine:
    input:
        lambda wildcards: expand("output/{run}/freebayes.vcf", run = SAMPLES[wildcards.sample])
    output:
        "output/{sample}/freebayes_combined.vcf"
    resources:
        runtime = 120,
        mem_mb = 4000
    wrapper:
        WRAPPER_PREFIX + "master/vcflib/vcfcombine"


rule vcffilter:
    input:
        "output/{sample}/freebayes_combined.vcf"
    output:
        "output/{sample}/freebayes_filtered.vcf"
    params:
        extra = "-f 'QUAL > 30 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1'"
    resources:
        runtime = 120,
        mem_mb = 4000
    wrapper:
        WRAPPER_PREFIX + "master/vcflib/vcffilter"


rule referencemaker:
    input:
        vcf = "output/{sample}/freebayes_combined.vcf",
        ref = REF_GENOME
    output:
        idx = temp("output/{sample}/freebayes_combined.vcf.idx"),
        fasta = "output/{sample}/consensus.fa",
        dic = "output/{sample}/consensus.dict",
        fai = "output/{sample}/consensus.fa.fai"
    params:
        refmaker = "--lenient"
    resources:
        runtime = 120,
        mem_mb = 4000    
    wrapper:
        WRAPPER_PREFIX + "master/gatk/fastaalternatereferencemaker"


rule fixfastaheader:
    input:
        rules.referencemaker.output.fasta
    output:
        "output/{sample}/consensus_fix.fa"
    params:
        run = lambda wildcards: wildcards.sample
    resources:
        runtime = 120,
        mem_mb = 2000
    shell:
        "sed 's/>.*/>{params.sample}/' {input[0]} > {output}"


rule sarscov2seqs:
    output:
        "output/sars-cov-2/sequences.gb"
    params:
        email = "taavi.pall@ut.ee",
        api_key = os.environ.get("NCBI_APIKEY")
    resources:
        runtime = 120,
        mem_mb = 2000
    wrapper:
        "file:wrappers/sequences/get_gb"


rule parsegb:
    input:
        "output/sars-cov-2/sequences.gb"
    output:
        fasta = "output/sars-cov-2/sequences.fa",
        metadata = "output/sars-cov-2/metadata.tsv"
    resources:
        runtime = 120,
        mem_mb = 2000
    wrapper:
        "file:wrappers/sequences/parse_gb"


# Run cd-hit to cluster identical sequences
# Shorter sequences will be clustered with longer ones if 
# they are completely covered by longer one
rule cd_hit:
    input:
        rules.parsegb.output.fasta
    output:
        repres = "output/sars-cov-2/cdhit.fa",
        clstr = "output/sars-cov-2/cdhit.fa.clstr"
    params:
        extra = "-c 1 -G 0 -aS 1 -d 0 -M 0"
    log:
        "output/sars-cov-2/log/cdhit.log"
    threads: 4
    resources:
        runtime = 120,
        mem_mb = 4000
    wrapper:
        WRAPPER_PREFIX + "master/cdhit"


rule align:
    input:
        "output/sars-cov-2/cdhit.fa"
    output:
        "output/sars-cov-2/msa.fa"
    log:
       "output/sars-cov-2/log/mafft.log"
    resources:
        runtime = 120,
        mem_mb = 4000
    wrapper:
        "file:wrappers/mafft"


rule merge:
    input:
        "output/{sample}/consensus_fix.fa",
        "output/sars-cov-2/msa.fa"
    output:
        "output/{sample}/msa.fa"
    params:
        extra="--add"
    log:
       "output/{sample}/log/mafft.log"
    resources:
        runtime = 120,
        mem_mb = 4000
    wrapper:
        "file:wrappers/mafft"
