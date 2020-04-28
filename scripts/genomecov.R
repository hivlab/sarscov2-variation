#+ include=FALSE
knitr::opts_chunk$set(echo = FALSE, message=FALSE)

#+ libs
library(tidyverse)
library(plotly)

#+ files
files <- list.files(here::here("output/genomecov"), recursive = TRUE, full.names = TRUE)

#+ map-files
cov <- map(files, read_tsv, col_names = c("chrom", "start", "end", "n"))

#+ runs
names(cov) <- str_split(files, "/") %>% map_chr(~ .x[length(.x) - 1])

#+ cov-long
cov_long <- cov %>% 
  bind_rows(.id = "Run") %>% 
  mutate(Sample = str_replace(Run, "-I", ""),
         Amplicons = case_when(
           str_detect(Run, "-I") ~ "set I",
           TRUE ~ "set non-I"
         )) %>% 
  select(Sample, Run, everything())

#+ plot
p <- cov_long %>% 
  pivot_longer(c("start",	"end"), values_to = "pos") %>% 
  ggplot() +
  geom_line(aes(pos, n, group = Run, color = Amplicons)) +
  facet_wrap(~Sample) +
  scale_y_log10()
ggplotly(p)

basecoverage <- cov_long %>% 
  mutate(pos = map2(start, end, ~seq(.x, .y - 1))) %>% 
  select(Sample, Run, pos, n) %>% 
  unnest(cols = pos)
basecoverage %>% 
  group_by(Sample, pos) %>% 
  summarise(n = sum(n)) %>% 
  filter(n > 10) %>% 
  group_by(Sample) %>% 
  count() %>% 
  arrange(desc(n)) %>% 
  mutate(`%` = (n / 29903) * 100)
