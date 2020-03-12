#' ---
#' title: "SARS-CoV-2 sequencing report"
#' author: "`r author`"  
#' date: "`r Sys.Date()`"
#' ---

#+ opts, include=FALSE
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE, fig.align='center')

#+ libs
library(tidyverse)
library(here)
library(knitr)
library(kableExtra)
library(formattable)
library(glue)
library(plotly)

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
  map(read_tsv, col_names = FALSE)

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
rl_plot <- ggplot(read_lengths) +
  geom_col(aes(X2, X3)) +
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
cov_plot <- cov_dist %>% 
  ggplot(aes(x=X3, y=X4)) +
  geom_point() +
  geom_line(aes(group = 1)) +
  labs(x = "Coverage", y = "Sites") +
  scale_x_log10()
ggplotly(cov_plot)

#' ## Variant calling
vcf <- read_tsv(vcf_file, comment = "#", col_names = FALSE)
colnames(vcf) <- str_split("CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	output/SRR11092064.bam", pattern = "\\s+", simplify = TRUE)
vcf %>% 
  mutate_at("POS", formatC, big.mark=",", format="d") %>% 
  select(CHROM:INFO) %>% 
  kable(caption = "Variant calling summary.") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = FALSE)




