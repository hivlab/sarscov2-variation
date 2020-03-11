
wd <- getwd()

# params
author <- snakemake@params[["author"]]
run <- snakemake@params[["run"]]

# iofiles
bamstats_file <- file.path(wd, snakemake@input[["bamstats"]])
vcf_file <- file.path(wd, snakemake@input[["vcf"]])
output_file <- file.path(wd, snakemake@output[[1]])

# render
rmarkdown::render(
  here::here("wrappers/report/report.R"), 
  output_file = output_file
  )
