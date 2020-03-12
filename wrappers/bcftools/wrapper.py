__author__ = "Taavi Päll"
__copyright__ = "Copyright 2020, Taavi Päll"
__email__ = "taavi.pall@ut.ee"
__license__ = "MIT"


from snakemake.shell import shell


shell(
    """
    bcftools mpileup {snakemake.params.mpileup} -f {snakemake.input.ref} {snakemake.input.samples} | \
    bcftools call {snakemake.params.call} | \
    bcftools filter {snakemake.params.filter} > {snakemake.output[0]}
    """
)
