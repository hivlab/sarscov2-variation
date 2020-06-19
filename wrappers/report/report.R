#' ---
#' title: "SARS-CoV-2 sequencing report"
#' author: "`r author`"  
#' date: "`r Sys.Date()`"
#' ---


#+ include=FALSE
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE, fig.align='center', knitr.duplicate.label = 'allow')


#+ libs
pkg <- c("dplyr", "readr", "purrr", "tidyr", "stringr", "here", "knitr", "kableExtra", "formattable", "glue", "plotly", "DT", "jcolors")
invisible(lapply(pkg, library, character.only = TRUE))


#+ options



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
statsfile <- read_lines(statsfile)
reads_used <- statsfile[1] %>% 
  str_split("\\t", simplify = TRUE) %>% 
  str_trim() %>% 
  str_c(collapse = " ") %>% 
  str_to_sentence() %>% 
  str_c(".")
summary_nums <- statsfile[str_which(statsfile, "mapped"):length(statsfile)] %>% 
  str_subset("^(?![\\s\\S])", negate = TRUE) %>% 
  str_c(collapse = "\n") %>% 
  read_tsv(col_names = FALSE)
colnames(summary_nums) <- c("variable", "% reads", "num reads", "% bases", "num bases")
summary_nums %>%
  kable(caption = glue("Summary of alignment of {run} to refseq. {reads_used}")) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = FALSE)


#' ### GC content
#+ gc-first-fragments, fig.cap='GC content of first fragments.'
gc_cont <- read_delim(gchist, comment = "#", delim = "\t", col_names = c("GC",	"Count"))
gc_plot <- gc_cont %>% 
  filter(Count > 0) %>% 
  mutate(Fraction = (Count / sum(Count)) * 100) %>% 
  ggplot() +
  geom_line(aes(GC, Fraction, group = 1)) +
  labs(x = "GC%", y = "Reads %")
ggplotly(gc_plot)


#' ### Average read quality
#+ aqhist, fig.cap='Histogram of average read quality.'
aq_reads <- read_delim(aqhist, comment = "#", delim = "\t", col_names = c("Quality", "Count", "Fraction"))
aq_plot <- aq_reads %>% 
  filter(Count > 0) %>% 
  mutate(Fraction = Fraction * 100) %>% 
  ggplot() +
  geom_line(aes(Quality, Fraction, group = 1)) +
  labs(x = "Average quality", y = "Reads %")
ggplotly(aq_plot)


#' ### Read lengths
#+ read-length, fig.cap='Read lengths.'
read_length <- read_delim(lhist, comment = "#", delim = "\t", col_names = c("Length", "Count"))
length_plot <- read_length %>% 
  filter(Count > 0) %>% 
  mutate(Fraction = (Count / sum(Count)) * 100) %>% 
  ggplot() +
  geom_line(aes(Length, Fraction, group = 1)) +
  labs(y = "Reads %")
ggplotly(length_plot)


#' ### Base composition
#+ base-comp, fig.cap='Distribution of base composition by position.'
base_comp <- read_delim(bhist, comment = "#", delim = "\t", col_names = c("Pos", "A", "C", "G", "T", "N"))
base_plot <- base_comp %>% 
  pivot_longer(c("A", "C", "G", "T", "N")) %>% 
  ggplot() +
  geom_line(aes(Pos, value, color = name, group = name)) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(y = "Reads %", x = "Base position") +
  theme(legend.title = element_blank(),
        legend.position = "bottom") +
  scale_color_jcolors(palette = "pal7")
base_plot


#' ### Mutations distribution
#+ mut-dist, fig.cap='Mutations distribution.'
indel_dist <- read_delim(mhist, comment = "#", delim = "\t", col_names = c("BaseNum",	"Match",	"Sub",	"Del",	"Ins",	"N",	"Other"))
indel_plot <- indel_dist %>% 
  pivot_longer(c("Match",	"Sub",	"Del",	"Ins",	"N",	"Other")) %>% 
  mutate_at("name", factor, levels = c("Match",	"Sub",	"Del",	"Ins",	"N",	"Other")) %>% 
  ggplot() +
  geom_line(aes(BaseNum, value, color = name, group = name)) +
  labs(y = "Reads %", x = "Base position") +
  theme(legend.title = element_blank(),
        legend.position = "bottom") +
  scale_color_jcolors(palette = "pal7") +
  scale_y_continuous(labels = scales::percent)
indel_plot

#' ### Coverage distribution
#+ cov-dist, fig.cap='Distribution of the alignment depth per covered reference site.'
genomecov <- read_delim(genomecov, comment = "#", delim = "\t", col_names = c("Ref",	"Start",	"End",	"Coverage"))
cov_plot <- genomecov %>% 
  pivot_longer(c("Start",	"End"), values_to = "Pos") %>% 
  ggplot(aes(Pos, Coverage)) +
  geom_line() +
  scale_y_log10() +
  scale_x_continuous(limits = c(0, 29903))
ggplotly(cov_plot)

#' ## Sequence variants
#+
parse_vcf_info <- function(x) {
  description <- map(x, str_extract, '(?<=Description=\\").+?(?=\\")')
  names(description) <- map(x, str_extract, "(?<=ID\\=)\\w+")
  str_c(str_c(names(description), str_remove(str_to_lower(description), "\\.$"), sep = ", ", collapse = ". "), ".")
}
# vcf_file <- "test/output/SRR11140750/freebayes.vcf"
vcf <- read_tsv(vcf_file, comment = "##")
comments <- read_lines(vcf_file) %>% 
  grep("##", ., value = TRUE)
info <- c("##FILTER", "##INFO", "##FORMAT") %>% 
  map(~grep(., comments, value = TRUE)) %>% 
  map(parse_vcf_info) %>% 
  str_c(collapse = " ")

parse_info <- function(data) {
  data <- str_split(data, ";")
  data = map(data, str_split, "=")
  data = map(data, transpose)
  data = map(data, ~map(.x, unlist))
  read_delim(
    str_c(
      str_c(data[[1]][[1]], collapse = " "), 
      str_c(data[[1]][[2]], collapse = " "), 
      sep = "\n"
      ), 
    delim = " ", 
    col_types = "dddc"
    )
}

vcf %>% 
  mutate(INFO = map(INFO, parse_info),
         ALT = if_else(str_length(ALT) > 10, str_c(str_sub(ALT, 1, 10), "..."), ALT)) %>%
  unnest(INFO) %>% 
  datatable(
    rownames = FALSE, 
    caption = str_c("Table 1. Variant calling summary. ", info),
    filter = "top", 
    options = list(
      pageLength = 20, 
      autoWidth = TRUE
    )
  )
