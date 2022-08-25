#!/usr/bin/env Rscript
# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
# Create a parser
p <- arg_parser("a program for QTLseq, TABEL file from GATK as input, 0 parent is needed")

## Add command line arguments
#
p <- add_argument(p, "--list", short = "-l", help = "list table header", flag = TRUE)
#
p <- add_argument(p, "--input", help = "TABEL file from GATK", type = "character")
p <- add_argument(p, "--out", help = "A prefix for output file", type = "character")
p <- add_argument(p, "--highB", help = "The bulk name with high phenotype", type = "character")
p <- add_argument(p, "--lowB", help = "The bulk name with low phenotype", type = "character")
#
p <- add_argument(p, "--minQ", help = "Minimum QUAL for variation", type = "numeric", default = 30)
#
#p <- add_argument(p, "--bulkSizeH", help = "Bulk size with high phenotype", type = "numeric")
#p <- add_argument(p, "--bulkSizeL", help = "Bulk size with low phenotype", type = "numeric")
#
p <- add_argument(p, "--minHBdp", help = "Minimum depth for high bulk", type = "numeric", default = 6)
p <- add_argument(p, "--maxHBdp", help = "Maxmum depth for high bulk", type = "numeric", default = 150)
p <- add_argument(p, "--minLBdp", help = "Minimum depth for low bulk", type = "numeric", default = 6)
p <- add_argument(p, "--maxLBdp", help = "Maxmum depth for low bulk", type = "numeric", default = 150)
#
p <- add_argument(p, "--winSize", help = "Window size for sliding window statistics", type = "numeric", default = 1000000)
p <- add_argument(p, "--winStep", help = "Window step for sliding window statistics", type = "numeric", default = 200000)

#
p <- add_argument(p, "--minN", help = "Minimum SNP number in a window for plot", type = "numeric", default = 10)
#
p <- add_argument(p, "--width", help = "Delta SNP index plot width", type = "numeric", default = 15)
p <- add_argument(p, "--height", help = "Delta SNP index plot height", type = "numeric", default = 5)

# Parse the command line arguments
argv <- parse_args(p)
#
if (argv$list) {
  df <- read.table(argv$input, nrows = 1, header = T)
  col_name <- colnames(df)
  print(col_name)
  quit(save = "no")
}
write.table(argv, paste(argv$out, "QTLseqrparameter.txt", sep = "."), quote = F, sep = "\t", row.names = F)

library(tidyverse)
#library(QTLseqr)
library(windowscanr)
library(cowplot)
#library(CMplot)
library(ggsci)
library(RColorBrewer)

filename <- argv$input
outPrefix <- argv$out
minQ <- argv$minQ
#bulkSizeH <- argv$bulkSizeH
#bulkSizeL <- argv$bulkSizeL
highB <- argv$highB
lowB <- argv$lowB
minHBdp <- argv$minHBdp
maxHBdp <- argv$maxHBdp
minLBdp <- argv$minLBdp
maxLBdp <- argv$maxLBdp
winSize <- argv$winSize
winStep <- argv$winStep
minN <- argv$minN
width <- argv$width
height <- argv$height

if (FALSE) {
  ## file name
  filename <- "./purpleLeaf.filter.SNPs.table.gz"
  outPrefix <- "purpleLeaf"
  ## filter SNP Quality
  minQ <- 30
  ##
  #bulkSizeH <- 30
  #bulkSizeL <- 30
  ## sample name
  highB <- "F2P"
  lowB <- "F2G"
  ##
  ## filter depth parameter
  minHBdp <- 6
  maxHBdp <- 150
  
  minLBdp <- 6
  maxLBdp <- 150
  
  ## sliding window parameter
  winSize <- 2000000
  winStep <- 200000
  #
  minN <- 10
  ##
  width <- 15
  height <- 5
}

## 读取数据
chr <- read_tsv("./chrom.txt")
len <- read_tsv(file = "./ref.len", col_names = c("CHROM", "Len"))
source("./Support_functions.R")

df <- read_tsv(file = filename) %>%
  select(CHROM, POS, REF, ALT, 
         starts_with(all_of(highB)), 
         starts_with(all_of(lowB))) %>% 
  na.omit() %>% rename(highBulk.GT = paste(highB, "GT", sep = "."),
                       lowBulk.GT = paste(lowB, "GT", sep = "."),
                       highBulk.AD = paste(highB, "AD", sep = "."),
                       lowBulk.AD = paste(lowB, "AD", sep = "."),
                       highBulk.DP = paste(highB, "DP", sep = "."),
                       lowBulk.DP = paste(lowB, "DP", sep = "."),
                       highBulk.GQ = paste(highB, "GQ", sep = "."),
                       lowBulk.GQ = paste(lowB, "GQ", sep = "."),
                       highBulk.PL = paste(highB, "PL", sep = "."),
                       lowBulk.PL = paste(lowB, "PL", sep = "."))

## 其实直接过滤两个混池的SNP index<0.3或>0.7的也行，只不过现在先过滤一部分有助于后面的运算速度提升
# 这里有潜在的危险，highBulk.AD_0 + highBulk.AD_1有时候小于highBulk.DP，这里过滤highBulk.DP后，
# 后面highBulk.AD_0 + highBulk.AD_1依然有潜在等于0的可能，导致后面做分母求index是报错（可能性很小，暂时不管了）
dd1 <- df %>% filter(!((highBulk.GT == paste(REF, REF, sep = "/") & lowBulk.GT == paste(REF, REF, sep = "/")) | 
                        (highBulk.GT == paste(ALT, ALT, sep = "/") & lowBulk.GT == paste(ALT, ALT, sep = "/")))) %>%
  filter(highBulk.DP > minHBdp, highBulk.DP < maxHBdp, lowBulk.DP > minLBdp, lowBulk.DP < maxLBdp)
dd2 <- dd1 %>% separate(highBulk.AD, c("highBulk.AD_0", "highBulk.AD_1"), sep = ",", convert = TRUE) %>%
  separate(lowBulk.AD, c("lowBulk.AD_0", "lowBulk.AD_1"), sep = ",", convert = TRUE) %>% 
  mutate(highBulk.ref.index = highBulk.AD_0 / (highBulk.AD_0 + highBulk.AD_1),
         highBulk.alt.index = highBulk.AD_1 / (highBulk.AD_0 + highBulk.AD_1),
         lowBulk.ref.index = lowBulk.AD_0 / (lowBulk.AD_0 + lowBulk.AD_1),
         lowBulk.alt.index = lowBulk.AD_1 / (lowBulk.AD_0 + lowBulk.AD_1),
         ED = sqrt((highBulk.ref.index - lowBulk.ref.index)^2 + (highBulk.alt.index - lowBulk.alt.index)^2))
dd3 <- dd2 %>% filter(!((highBulk.ref.index < 0.3 & lowBulk.ref.index < 0.3) |
                          (highBulk.ref.index > 0.7 & lowBulk.ref.index > 0.7)))
## 
dp <- dd3 %>% dplyr::select(HB = highBulk.DP, LB = lowBulk.DP) %>%
  gather(key = "sample", value = "depth")
dp$sample <- factor(dp$sample, levels = c("HB", "LB"), labels = c(highB, lowB))
P_dp <- ggplot(dp, aes(x = depth)) + 
  geom_histogram(aes(y = ..density.., fill = sample), binwidth = 2) +
  geom_density() + 
  scale_x_continuous(limits = c(0, 100)) +
  theme_half_open() +
  facet_wrap(~sample, nrow = 1)
ggsave(P_dp, filename = paste(outPrefix, "depth_desity.pdf", sep = "."), height = 3.5, width = 6)
ggsave(P_dp, filename = paste(outPrefix, "depth_desity.png", sep = "."), height = 3.5, width = 6, dpi = 500)
remove(P_dp)

###
##
SNPnumber <- dd3 %>% group_by(CHROM) %>% count()
write_tsv(x = SNPnumber, file = paste(outPrefix, "SNP_number_per_chr.txt", sep = "."))
write_csv(x = SNPnumber, file = paste(outPrefix, "SNP_number_per_chr.csv", sep = "."))
##
options(scipen = 200)
colourCount = dim(chr)[[1]]
getPalette = colorRampPalette(brewer.pal(8, "Set1"))
#df$LABEL <- factor(df$LABEL, levels = chromColor$LABEL)
chr$LABEL <- factor(chr$LABEL, levels = chr$LABEL)
Phist <- chr %>% left_join(dd4, by = "CHROM") %>% ggplot(aes(x = POS)) +
  geom_histogram(aes(fill = LABEL), color = NA, binwidth = 1000000) +
  labs(x = NULL, y = "SNP Count / 1Mb", fill = "Chrom") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(n.breaks = 2.1) +
  scale_fill_manual(values = getPalette(colourCount)) +
  theme_half_open() +
  theme(strip.text.y = element_text(angle = 0),
        strip.background = element_rect(color = NA, fill = NA),
        legend.position = "NULL") +
  facet_grid(LABEL ~ .)
ggsave(Phist, filename = paste(outPrefix, "SNP_distribution_histogram.pdf", sep = "."), width = 8, height = dim(chr)[[1]] * 0.45 + 0.5)
ggsave(Phist, filename = paste(outPrefix, "SNP_distribution_histogram.png", sep = "."), width = 8, height = dim(chr)[[1]] * 0.45 + 0.5, dpi = 500)
remove(Phist)

##
slidwin <- winScan(x = dd3,
                   groups = "CHROM",
                   position = "POS",
                   values = c("ED"),
                   win_size = winSize,
                   win_step = winStep,
                   funs = c("mean")) %>%
  dplyr::select(CHROM, win_start, win_end, win_mid, nSNPs = ED_n, ED = ED_mean) %>%
  dplyr::mutate(ED4 = ED^4)
write_tsv(x = slidwin, file = paste(outPrefix, "SlidingWindow.txt", sep = "."))
write_csv(x = slidwin, file = paste(outPrefix, "SlidingWindow.csv", sep = "."))

slidwin$POS <- slidwin$win_mid
## ED
pdf(file = paste(outPrefix, "ED.point.pdf", sep = "."), width = width, height = height)
plotIndex(df = slidwin, chr = chr, len = len, nSNPs = "nSNPs", n = minN, 
          index = "ED", band = 0.005, ylim = c(0, max(slidwin$ED, na.rm = TRUE)),
          ylab = "ED", size = 1, axis.size = 1, 
          axis.lab.size = 1, axis.title.size = 1.5, type = "p")
dev.off()
png(filename = paste(outPrefix, "ED.point.png", sep = "."), width = width, height = height, units = "in", res = 500)
plotIndex(df = slidwin, chr = chr, len = len, nSNPs = "nSNPs", n = minN, 
          index = "ED", band = 0.005, ylim = c(0, max(slidwin$ED, na.rm = TRUE)),
          ylab = "ED", size = 1, axis.size = 1, 
          axis.lab.size = 1, axis.title.size = 1.5, type = "p")
dev.off()
## ED^4
pdf(file = paste(outPrefix, "ED4.point.pdf", sep = "."), width = width, height = height)
plotIndex(df = slidwin, chr = chr, len = len, nSNPs = "nSNPs", n = minN, 
          index = "ED4", band = 0.005, ylim = c(0, max(slidwin$ED4, na.rm = TRUE)),
          ylab = "ED^4", size = 1, axis.size = 1, 
          axis.lab.size = 1, axis.title.size = 1.5, type = "p")
dev.off()
## ED
pdf(file = paste(outPrefix, "ED.line.pdf", sep = "."), width = width, height = height)
plotIndex(df = slidwin, chr = chr, len = len, nSNPs = "nSNPs", n = minN, 
          index = "ED", band = 0.005, ylim = c(0, max(slidwin$ED, na.rm = TRUE)),
          ylab = "ED", size = 1.5, axis.size = 1, 
          axis.lab.size = 1, axis.title.size = 1.5, type = "l")
dev.off()
png(filename = paste(outPrefix, "ED.line.png", sep = "."), width = width, height = height, units = "in", res = 500)
plotIndex(df = slidwin, chr = chr, len = len, nSNPs = "nSNPs", n = minN, 
          index = "ED", band = 0.005, ylim = c(0, max(slidwin$ED, na.rm = TRUE)),
          ylab = "ED", size = 1.5, axis.size = 1, 
          axis.lab.size = 1, axis.title.size = 1.5, type = "l")
dev.off()
## ED^4
pdf(file = paste(outPrefix, "ED4.line.pdf", sep = "."), width = width, height = height)
plotIndex(df = slidwin, chr = chr, len = len, nSNPs = "nSNPs", n = minN, 
          index = "ED4", band = 0.005, ylim = c(0, max(slidwin$ED4, na.rm = TRUE)),
          ylab = "ED^4", size = 1.5, axis.size = 1, 
          axis.lab.size = 1, axis.title.size = 1.5, type = "l")
dev.off()
png(filename = paste(outPrefix, "ED4.line.png", sep = "."), width = width, height = height, units = "in", res = 500)
plotIndex(df = slidwin, chr = chr, len = len, nSNPs = "nSNPs", n = minN, 
          index = "ED4", band = 0.005, ylim = c(0, max(slidwin$ED4, na.rm = TRUE)),
          ylab = "ED^4", size = 1.5, axis.size = 1, 
          axis.lab.size = 1, axis.title.size = 1.5, type = "l")
dev.off()

