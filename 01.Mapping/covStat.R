#!/bin/env Rscript
# need genome.fa.fai as argv
library(tidyverse)

argv <- commandArgs(trailingOnly = TRUE)
chr <- read_tsv(file = argv[[1]], col_names = FALSE)
genoSize <- chr %>% pull(X2) %>% sum()

all <- tibble()
#sample <- "test1"
for (sample in list.files(pattern = ".0cov.bedgraph") %>% basename() %>% str_remove(".0cov.bedgraph")) {
  cat(date(), ": reading ", sample, ".0cov.bedgraph ......", sep = "")
  df <- read_tsv(paste(sample, "0cov.bedgraph", sep = "."), col_names = FALSE)
  stat <- df %>% mutate(len = X3 - X2) %>% 
    summarise(covLen = genoSize - sum(len)) %>% 
    mutate(sample = sample, covRate = covLen / genoSize) %>%
    select(sample, covLen, covRate)
  all <- rbind(all, stat)
}

write_csv(x = all, path = "cov_stat.csv")
write_tsv(x = all, path = "cov_stat.txt")
