__author__ = "Taavi Päll"
__copyright__ = "Copyright 2020, Taavi Päll"
__email__ = "taavi.pall@ut.ee"
__license__ = "MIT"


from snakemake.shell import shell

# Parse params
extra = snakemake.params.get("extra", "")

# Parse inputs
ind = snakemake.input
inputs = " ".join(list({"-" + k + " " + v for k, v in ind.items()}))


shell("bedtools genomecov {extra} {inputs} > {snakemake.output[0]}")
