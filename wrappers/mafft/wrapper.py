__author__ = "Taavi Päll"
__copyright__ = "Copyright 2020, Taavi Päll"
__email__ = "tapa741@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell

# Extra arguments
extra = snakemake.params.get("extra", "")

# Setup logfile
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# Run command
shell(
    "(mafft {extra} {snakemake.input} > {snakemake.output[0]}) {log}"
)
