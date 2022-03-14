#!/usr/bin/env Rscript
# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
# Create a parser
p <- arg_parser("extract functional annotation")

## Add command line arguments
#
p <- add_argument(p, "--bed", help = "bed file name", type = "character")
p <- add_argument(p, "--ann", help = "annotation file name", type = "character")
p <- add_argument(p, "--chr", help = "chromosome", type = "character")
p <- add_argument(p, "--start", help = "interval start position", type = "numeric")
p <- add_argument(p, "--end", help = "interval end position", type = "numeric")

# Parse the command line arguments
argv <- parse_args(p)

bedfile <- argv$bed
annfile <- argv$ann
CHR <- argv$chr
START <- argv$start
END <- argv$end


test <- F
if (test) {
  bedfile <- "no2127.v0.bed"
  annfile <- "no2127.geneinformation.txt"
  CHR <- "no2127v0_C03"
  START <- "0"
  END <- "12446968"
}

library(tidyverse)

bed <- read_tsv(file = bedfile, col_names = F) %>% 
  select(chr = X1, start = X2, end = X3, Locus_ID = X4)
bed
ann <- read_tsv(file = annfile)
ann
df <- bed %>% left_join(ann, by = "Locus_ID") %>% 
  filter(chr == CHR & start < END & end > START)
outprefix <- paste(paste(CHR, START, END, sep = "_"))
write_csv(x = df, file = paste(outprefix, "csv", sep = "."))
write_tsv(x = df, file = paste(outprefix, "txt", sep = "."))
