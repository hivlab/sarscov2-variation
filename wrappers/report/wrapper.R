
rmarkdown::render(
  here::here("wrappers/report/report.R"), 
  output_file = here::here(snakemake@output[[1]]),
  params = list(run = snakemake@params[["run"]], 
                bamstats = snakemake@input[["bamstats"]],
                vcf = snakemake@input[["vcf"]],
                author = snakemake@params[["author"]])
  )
