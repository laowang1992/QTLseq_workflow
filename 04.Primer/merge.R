library(tidyverse)
library(xlsx)

argv <- commandArgs(trailingOnly = T)

gt <- read_tsv(file = argv[1])
primer <- read_tsv(file = "xxx.primer.txt")

df <- gt %>% right_join(primer, by = "ID")
write.csv(x = df, file = str_replace(argv[1], "primer.txt", "primer.csv"))
write.xlsx(x = df, file = str_replace(argv[1], "primer.txt", "primer.xlsx"))
