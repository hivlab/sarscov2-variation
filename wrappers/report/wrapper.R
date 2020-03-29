
wd <- getwd()

# params
author <- snakemake@params[["author"]]
run <- snakemake@params[["run"]]

# iofiles
statsfile <- file.path(wd, snakemake@input[["statsfile"]])
gchist <- file.path(wd, snakemake@input[["gchist"]])
aqhist <- file.path(wd, snakemake@input[["aqhist"]])
lhist <- file.path(wd, snakemake@input[["lhist"]])
mhist <- file.path(wd, snakemake@input[["mhist"]])
bhist <- file.path(wd, snakemake@input[["bhist"]])
genomecov <- file.path(wd, snakemake@input[["genomecov"]])
vcf_file <- file.path(wd, snakemake@input[["vcf"]])
output_file <- file.path(wd, snakemake@output[[1]])

# render
rmarkdown::render(
  here::here("wrappers/report/report.R"), 
  output_file = output_file,
  output_format = bookdown::html_document2(number_sections = FALSE, theme = "flatly"),
  envir = new.env()
  )
