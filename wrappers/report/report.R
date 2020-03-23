#' ---
#' title: "SARS-CoV-2 sequencing report"
#' author: "`r author`"  
#' date: "`r Sys.Date()`"
#' ---

#+ include=FALSE
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE, fig.align='center')

#+ libs
pkg <- c("tidyverse", "here", "knitr", "kableExtra", "formattable", "glue", "plotly", "DT")
invisible(lapply(pkg, library, character.only = TRUE))

#+ parse
bamstats <- read_lines(bamstats_file)
vars <- str_subset(bamstats, "^#.*Use") %>%
  str_extract("^[^.]+") %>%
  str_replace("# ", "")
var_regex <- bamstats %>%
  str_subset("^#.*Use") %>%
  str_extract("\\^[A-Z]{2,3}")
values <- var_regex %>%
  map(~str_subset(bamstats, .x)) %>%
  map(str_c, collapse="\n")
names(values) <- vars %>%
  str_to_lower() %>%
  str_replace_all("[ -]+", "_")
parsed_stats <- values %>%
  keep(~str_length(.x) != 0) %>%
  map(~as_tibble(read.csv(text = .x, header = FALSE, sep = "\t")))

#' ## Refseq
#' 
#' - NCBI Reference Sequence: *NC_045512.2* -- Wuhan seafood market pneumonia virus isolate Wuhan-Hu-1, complete genome
#'    - A novel coronavirus associated with a respiratory disease in Wuhan of Hubei province, China
#'    - 29903 bp, ss-RNA, linear, VRL, 28-JAN-2020
#'

#+ bam-stats-chapter, results='asis'
tabset <- "{.tabset .tabset-fade}"
glue("## {run} alignment to refseq {tabset}")
#' 
 

#' ### Summary numbers
#+ summary-nums
summary_nums <- parsed_stats[["summary_numbers"]] %>% select(-1)
colnames(summary_nums) <- c("key", "value")
summary_nums %>%
  filter(value > 0) %>% 
  mutate_at("value", formatC, big.mark=",", format="d") %>% 
  kable(caption = "Summary numbers.") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = FALSE)

#' ### GC content
#+ gc-first-fragments, fig.cap='GC content of first fragments.'
gc_cont <- parsed_stats[["gc_content_of_first_fragments"]]
gc_plot <- tibble(x=gc_cont[[2]], y=gc_cont[[3]]) %>% 
  ggplot() +
  geom_line(aes(x, y, group = 1)) +
  labs(x = "GC%", y = "First fragments")
ggplotly(gc_plot)


#' ### AGCT content
#+ acgt-content, fig.cap='AGCT content per cycle.'
agct_cycle <- parsed_stats[["acgt_content_per_cycle"]] %>% select(-1)
colnames(agct_cycle) <- c("cycle", LETTERS[c(1, 3, 7, 20)], "N", "O")
agct_plot <- agct_cycle %>% 
  pivot_longer(c("A", "C", "G", "T", "N", "O")) %>% 
  filter(name %in% c("A", "C", "G", "T")) %>% 
  ggplot() +
  geom_tile(aes(name, cycle, fill=value)) +
  labs(y = "Cycle") +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank())
ggplotly(agct_plot)

#' ### Read lengths
#+ read-length, fig.cap='Read lengths.'
read_lengths <- parsed_stats[["read_lengths"]] %>% select(-1)
colnames(read_lengths) <- c("read_length", "count")
rl_plot <- ggplot(read_lengths) +
  geom_col(aes(read_length, count)) +
  labs(x = "Read length", y = "Count") +
  scale_y_log10()
ggplotly(rl_plot)

#' ### In/del distribution
#+ indel-dist, fig.cap='In/del distribution.'
indel_dist <- parsed_stats[["indel_distribution"]] %>% select(-1)
colnames(indel_dist) <- c("length", "ins", "del")
indel_plot <- indel_dist %>% 
  pivot_longer(c("ins", "del")) %>% 
  mutate_at("name", factor, levels = c("ins", "del")) %>% 
  ggplot() +
  geom_col(aes(length, value)) +
  facet_wrap(~name) +
  labs(x = "Length", y = "Count")
ggplotly(indel_plot)

#' ### Coverage distribution
#+ cov-dist, fig.cap='Distribution of the alignment depth per covered reference site.'
cov_dist <- parsed_stats[["coverage_distribution"]] %>% select(-1)
colnames(cov_dist) <- c("range", "cov", "sites")
cov_plot <- cov_dist %>% 
  ggplot(aes(cov, sites)) +
  geom_point() +
  geom_line(aes(group = 1)) +
  labs(x = "Coverage", y = "Sites") +
  scale_x_log10()
ggplotly(cov_plot)

#' ## Sequence variants
#+
parse_vcf_info <- function(x) {
  description <- map(x, str_extract, '(?<=Description=\\").+?(?=\\")')
  names(description) <- map(x, str_extract, "(?<=ID\\=)\\w+")
  str_c(str_c(names(description), str_remove(str_to_lower(description), "\\.$"), sep = ", ", collapse = ". "), ".")
}
vcf <- read_tsv(vcf_file, comment = "##")
comments <- read_lines(vcf_file) %>% 
  grep("##", ., value = TRUE)
info <- c("##FILTER", "##INFO", "##FORMAT") %>% 
  map(~grep(., comments, value = TRUE)) %>% 
  map(parse_vcf_info) %>% 
  str_c(collapse = " ")
vcf %>% 
  select(`#CHROM`:INFO) %>% 
  datatable(
    rownames = FALSE, 
    caption = str_c("Table 1. Variant calling summary. ", info),
    filter = "top", 
    options = list(
      pageLength = 20, 
      autoWidth = TRUE
      )
    )
