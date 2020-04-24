library(tidyverse)

snps <- read_csv("output/snpsift.csv")

split_var <- function(x) {
  unlist(str_split(x, ","))
}

get_length <- function(x) {
  length(split_var(x)) > 1
}

map_split_var <- function(x) {
  map(x, split_var)
}

chunks <- function(x, n) {
  assertthat::assert_that((length(x) %% n) == 0, msg = "Length of input must be even.")
  if (length(x) > n) {
    d <- rep(1:n, each = length(x) / n)
    return(map_chr(split(x, d), paste0, collapse = ","))
  } 
  return(x)
}


multi <- snps %>%
  rename(Run = Sample) %>% 
  mutate(Sample = str_replace(Run, "-I", ""),
         multi_snps = map_lgl(ALT, get_length)) %>% 
  split(.$multi_snps)
one <- multi$"FALSE" %>% 
  mutate_at(vars(QUAL:RPL, AF), as.numeric)
more <- multi$"TRUE"
unnested <- more %>% 
  mutate_at(vars(ALT, QUAL, AO, SAF, SAR, RPR, RPL, AF, starts_with("ANN")), map_split_var) %>% 
  mutate(len_alt = map_int(ALT, length), 
         `ANN[*].IMPACT` = map2(`ANN[*].IMPACT`, len_alt, chunks), 
         `ANN[*].EFFECT` = map2(`ANN[*].EFFECT`, len_alt, chunks), 
         `ANN[*].GENE` = map2(`ANN[*].GENE`, len_alt, chunks), 
         `ANN[*].CODON` = map2(`ANN[*].CODON`, len_alt, chunks)) %>% 
  unnest(cols = c(ALT, QUAL, AO, SAF, SAR, RPR, RPL, AF, `ANN[*].IMPACT`, `ANN[*].EFFECT`, `ANN[*].GENE`, `ANN[*].CODON`)) %>% 
  mutate_at(vars(QUAL:RPL, AF), as.numeric) %>% 
  select(-multi_snps, -len_alt)
multi_fix <- bind_rows(one, unnested)
multi_fix %>% 
  select(Sample, Run, everything()) %>% 
  select(-level_1, -multi_snps) %>% 
  write_csv("output/snpsift_unnested.csv")

# - QUAL > 1: removes very bad sites
# - QUAL / AO > 10: additional contribution of each obs should be 10 log units (~ Q10 per read) 
# - SAF > 0 & SAR > 0: reads on both strands 
# - RPR > 1 & RPL > 1: at least two reads "balanced" to each side of the site


good_snps <- multi_fix %>% 
  select(Sample, Run, everything()) %>% 
  filter(QUAL > 30, QUAL / AO > 10, SAF > 0, SAR > 0) %>% 
  arrange(Sample, POS)

parsed <- good_snps %>% 
  select(Sample, Run, POS, ALT) %>% 
  distinct() %>% 
  group_by(Sample) %>% 
  summarise_at(vars(POS, ALT), str_c, collapse = ",")
parsed$POS[1]
