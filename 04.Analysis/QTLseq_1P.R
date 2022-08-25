#!/usr/bin/env Rscript
# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
# Create a parser
p <- arg_parser("a program for QTLseq, TABEL file from GATK as input, only one paren is needed")

## Add command line arguments
#
p <- add_argument(p, "--list", short = "-l", help = "list table header", flag = TRUE)
#
p <- add_argument(p, "--input", help = "TABEL file from GATK", type = "character")
p <- add_argument(p, "--out", help = "A prefix for output file", type = "character")

p <- add_argument(p, "--CI", help = "QTLseqCI**.RData, you must check whether it meet with you sample, e.g.: popType, bulk size, depth", type = "character", default = NA)

p <- add_argument(p, "--highP", help = "The parent name with high phenotype", type = "character")
p <- add_argument(p, "--lowP", help = "The parent name with low phenotype", type = "character")
p <- add_argument(p, "--highB", help = "The bulk name with high phenotype", type = "character")
p <- add_argument(p, "--lowB", help = "The bulk name with low phenotype", type = "character")
p <- add_argument(p, "--popType", help = "Population type, 'F2' or 'RIL'", type = "character")
#
p <- add_argument(p, "--minQ", help = "Minimum QUAL for variation", type = "numeric", default = 30)
#
p <- add_argument(p, "--bulkSizeH", help = "Bulk size with high phenotype", type = "numeric")
p <- add_argument(p, "--bulkSizeL", help = "Bulk size with low phenotype", type = "numeric")

p <- add_argument(p, "--minGQ", help = "Minimum GQ for parent genotype", type = "numeric", default = 30)
#
p <- add_argument(p, "--minHPdp", help = "Minimum depth for high parent", type = "numeric", default = 6)
p <- add_argument(p, "--maxHPdp", help = "Maxmum depth for high parent", type = "numeric", default = 150)
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

if (argv$list) {
  df <- read.table(argv$input, nrows = 1, header = T)
  col_name <- colnames(df)
  print(col_name)
  quit(save = "no")
}

write.table(argv, paste(argv$out, "QTLseqrparameter.txt", sep = "."), quote = F, sep = "\t", row.names = F)

#########################################
## load QTLseqCI data automatic
## not finish now, need load manually
if (FALSE) {
  minDepth <- min(minHBdp, minLBdp)
  maxDepth <- max(maxHBdp, maxLBdp)
  
  prefix <- paste("QTLseqCI", paste(paste("indH", bulkSizeH, sep = ""), paste("indL", bulkSizeL, sep = ""), popType, sep = "_"), sep = ".")
  CI_file <- list.files(pattern = paste("^", prefix, sep = ""))
  CI_info <- strsplit(x = CI, split = ".", fixed = TRUE)
}

library(tidyverse)
library(QTLseqr)
library(windowscanr)
library(cowplot)
#library(CMplot)
library(ggsci)
library(RColorBrewer)

filename <- argv$input
outPrefix <- argv$out
CI <- argv$CI
minQ <- argv$minQ
bulkSizeH <- argv$bulkSizeH
bulkSizeL <- argv$bulkSizeL
highP <- argv$highP
lowP <- argv$lowP
highB <- argv$highB
lowB <- argv$lowB
popType <- argv$popType
minGQ <- argv$minGQ
minHPdp <- argv$minHPdp
maxHPdp <- argv$maxHPdp
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
  #
  CI <- "./QTLseqCI.indH30_indL30_F2.Depth_5_200.Rep_10000.RData"
  ## filter SNP Quality
  minQ <- 30
  ##
  bulkSizeH <- 30
  bulkSizeL <- 30
  ## sample name
  highP <- "PP"
  lowP <- "PG"
  highB <- "F2P"
  lowB <- "F2G"
  ##
  popType <- "F2"  # F2 ro RIL
  ##
  minGQ <- 30
  ## filter depth parameter
  minHPdp <- 6
  maxHPdp <- 150
    
  minHBdp <- 6
  maxHBdp <- 150
  
  minLBdp <- 6
  maxLBdp <- 150
  
  ## sliding window parameter
  winSize <- 1000000
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
         starts_with(all_of(highP)), 
         starts_with(all_of(highB)), 
         starts_with(all_of(lowB))) %>% 
  na.omit() %>% rename(highParent.GT = paste(highP, "GT", sep = "."),
                       highBulk.GT = paste(highB, "GT", sep = "."),
                       lowBulk.GT = paste(lowB, "GT", sep = "."),
                       highParent.AD = paste(highP, "AD", sep = "."),
                       highBulk.AD = paste(highB, "AD", sep = "."),
                       lowBulk.AD = paste(lowB, "AD", sep = "."),
                       highParent.DP = paste(highP, "DP", sep = "."),
                       highBulk.DP = paste(highB, "DP", sep = "."),
                       lowBulk.DP = paste(lowB, "DP", sep = "."),
                       highParent.GQ = paste(highP, "GQ", sep = "."),
                       highBulk.GQ = paste(highB, "GQ", sep = "."),
                       lowBulk.GQ = paste(lowB, "GQ", sep = "."),
                       highParent.PL = paste(highP, "PL", sep = "."),
                       highBulk.PL = paste(highB, "PL", sep = "."),
                       lowBulk.PL = paste(lowB, "PL", sep = "."))

dd1 <- df %>% filter(highParent.GQ >= minGQ) %>% 
  filter(highParent.GT == paste(REF, REF, sep = "/") | highParent.GT == paste(ALT, ALT, sep = "/")) %>%
  filter(!((highBulk.GT == paste(REF, REF, sep = "/") | highBulk.GT == paste(ALT, ALT, sep = "/")) & 
             (highBulk.GT == lowBulk.GT)))

dd2 <- dd1 %>% separate(highBulk.AD, c("highBulk.AD_0", "highBulk.AD_1"), sep = ",", convert = TRUE) %>%
  separate(lowBulk.AD, c("lowBulk.AD_0", "lowBulk.AD_1"), sep = ",", convert = TRUE) %>%
  separate(highBulk.PL, c("highBulk.PL_00", "highBulk.PL_01", "highBulk.PL_11"), sep = ",", convert = TRUE) %>%
  separate(lowBulk.PL, c("lowBulk.PL_00", "lowBulk.PL_01", "lowBulk.PL_11"), sep = ",", convert = TRUE) %>%
  mutate(highBulk.DP = highBulk.AD_0 + highBulk.AD_1, lowBulk.DP = lowBulk.AD_0 + lowBulk.AD_1) # 因为highBulk.AD_0 + highBulk.AD_1有时并不等于原来的highBulk.DP并且会导致后面代码报错而且很难发现原因，因此这里重新计算覆盖原来的highBulk.DP

dd3 <- dd2 %>% filter(!(((highBulk.AD_0 / (highBulk.AD_0 + highBulk.AD_1) < 0.3) & 
                          (lowBulk.AD_0 / (lowBulk.AD_0 + lowBulk.AD_1) < 0.3)) | 
                        ((highBulk.AD_0 / (highBulk.AD_0 + highBulk.AD_1) > 0.7) & 
                           (lowBulk.AD_0 / (lowBulk.AD_0 + lowBulk.AD_1) > 0.7))))
## 
dp <- dd3 %>% dplyr::select(HP = highParent.DP, HB = highBulk.DP, LB = lowBulk.DP) %>%
  gather(key = "sample", value = "depth")
dp$sample <- factor(dp$sample, levels = c("HP", "HB", "LB"), labels = c(highP, highB, lowB))
P_dp <- ggplot(dp, aes(x = depth)) + 
  geom_histogram(aes(y = ..density.., fill = sample), binwidth = 2) +
  geom_density() + 
  scale_x_continuous(limits = c(0, 100)) +
  theme_half_open() +
  facet_wrap(~sample, nrow = 1)
ggsave(P_dp, filename = paste(outPrefix, "depth_desity.pdf", sep = "."), height = 3.5, width = 10)
ggsave(P_dp, filename = paste(outPrefix, "depth_desity.png", sep = "."), height = 3.5, width = 10, dpi = 500)
remove(P_dp)

dd4 <- dd3 %>% filter(highParent.DP > minHPdp, highParent.DP < maxHPdp,
         highBulk.DP > minHBdp, highBulk.DP < maxHBdp,
         lowBulk.DP > minLBdp, lowBulk.DP < maxLBdp)


##
SNPnumber <- dd4 %>% group_by(CHROM) %>% count()
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

## export QTLseqr input data
#dd5 <- dd4 %>% mutate(highBulk.AD = if_else(highParent.GT == paste(ALT, ALT, sep = "/"), 
#                                             paste(highBulk.AD_0, highBulk.AD_1, sep = ","), 
#                                             paste(highBulk.AD_1, highBulk.AD_0, sep = ",")),
#                       lowBulk.AD = if_else(highParent.GT == paste(ALT, ALT, sep = "/"),
#                                            paste(lowBulk.AD_0, lowBulk.AD_1, sep = ","),
#                                            paste(lowBulk.AD_1, lowBulk.AD_0, sep = ",")),
#                       highBulk.PL = if_else(highParent.GT == paste(ALT, ALT, sep = "/"),
#                                             paste(highBulk.PL_00, highBulk.PL_01, highBulk.PL_11, sep = ","),
#                                             paste(highBulk.PL_11, highBulk.PL_01, highBulk.PL_00, sep = ",")),
#                       lowBulk.PL = if_else(highParent.GT == paste(ALT, ALT, sep = "/"),
#                                            paste(lowBulk.PL_00, lowBulk.PL_01, lowBulk.PL_11, sep = ","),
#                                            paste(lowBulk.PL_11, lowBulk.PL_01, lowBulk.PL_00, sep = ","))) %>%
#  select(CHROM, POS, REF, ALT, 
#         highBulk.AD, highBulk.DP, highBulk.GQ, highBulk.PL,
#         lowBulk.AD, lowBulk.DP, lowBulk.GQ, lowBulk.PL)
#write_tsv(dd5, paste(outPrefix, "reform.txt", sep = "."))

dd6 <- dd4 %>% mutate(HB.HP.AD = if_else(highParent.GT == paste(REF, REF, sep = "/"),
                                         highBulk.AD_0, highBulk.AD_1),
                      HB.LP.AD = if_else(highParent.GT == paste(REF, REF, sep = "/"),
                                         highBulk.AD_1, highBulk.AD_0),
                      LB.HP.AD = if_else(highParent.GT == paste(REF, REF, sep = "/"),
                                         lowBulk.AD_0, lowBulk.AD_1),
                      LB.LP.AD = if_else(highParent.GT == paste(REF, REF, sep = "/"),
                                         lowBulk.AD_1, lowBulk.AD_0)) %>%
  dplyr::select(CHROM, POS, HB.HP.AD, HB.LP.AD, LB.HP.AD, LB.LP.AD) %>%
  dplyr::mutate(HB.DP = HB.HP.AD + HB.LP.AD, 
                LB.DP = LB.HP.AD + LB.LP.AD, 
                HB.index = HB.HP.AD / HB.DP, 
                LB.index = LB.HP.AD / LB.DP,
                delta.index = HB.index- LB.index,
                ED = sqrt((HB.index - LB.index)^2 + ((1-HB.index) - (1-LB.index))^2))

if (is.na(CI)) {
  cat("We will not give confidence interval information\n")
  slidwin <- winScan(x = dd6,
                     groups = "CHROM",
                     position = "POS",
                     values = c("HB.index", "LB.index", "delta.index", "ED"),
                     win_size = winSize,
                     win_step = winStep,
                     funs = c("mean")) %>% as_tibble() %>%
    dplyr::select(CHROM, win_start, win_end, win_mid, nSNPs = delta.index_n,
                  HB.index = HB.index_mean, LB.index = LB.index_mean, delta.index = delta.index_mean, 
                  ED = ED_mean) %>%
    dplyr::mutate(ED4 = ED^4)
  write_tsv(x = slidwin, file = paste(outPrefix, "SlidingWindow.txt", sep = "."))
  write_csv(x = slidwin, file = paste(outPrefix, "SlidingWindow.csv", sep = "."))
  
  ## 我们开始画图了
  #df <- chr %>% left_join(slidwin, by = "CHROM")
  #slidwin$CHROM <- factor(slidwin$CHROM, levels = chr$CHROM)
  ## SNP index
  # delta index
  slidwin$POS = slidwin$win_mid
  ####point plot
  ## Delta SNP index
  pdf(file = paste(outPrefix, "delta_SNP_index.point.pdf", sep = "."), width = width, height = height)
  plotIndex(df = slidwin, chr = chr, len = len, nSNPs = "nSNPs", n = minN, 
            index = "delta.index", band = 0.005, ylim = c(-1, 1),
            ylab = "Delta SNP index", size = 1, axis.size = 1, 
            axis.lab.size = 1, axis.title.size = 1.5, type = "p")
  dev.off()
  png(filename = paste(outPrefix, "delta_SNP_index.point.png", sep = "."), width = width, height = height, units = "in", res = 500)
  plotIndex(df = slidwin, chr = chr, len = len, nSNPs = "nSNPs", n = minN, 
            index = "delta.index", band = 0.005, ylim = c(-1, 1),
            ylab = "Delta SNP index", size = 1, axis.size = 1, 
            axis.lab.size = 1, axis.title.size = 1.5, type = "p")
  dev.off()
  ## High Bulk SNP index
  pdf(file = paste(outPrefix, highB, "SNP_index.point.pdf", sep = "."), width = width, height = height)
  plotIndex(df = slidwin, chr = chr, len = len, nSNPs = "nSNPs", n = minN, 
            index = "HB.index", band = 0.005, ylim = c(0, 1),
            ylab = paste(highB, "SNP index", sep = " "), size = 1, axis.size = 1, 
            axis.lab.size = 1, axis.title.size = 1.5, type = "p")
  dev.off()
  png(filename = paste(outPrefix, highB, "SNP_index.point.png", sep = "."), width = width, height = height, units = "in", res = 500)
  plotIndex(df = slidwin, chr = chr, len = len, nSNPs = "nSNPs", n = minN, 
            index = "HB.index", band = 0.005, ylim = c(0, 1),
            ylab = paste(highB, "SNP index", sep = " "), size = 1, axis.size = 1, 
            axis.lab.size = 1, axis.title.size = 1.5, type = "p")
  dev.off()
  ## Low Bulk SNP index
  pdf(file = paste(outPrefix, lowB, "SNP_index.point.pdf", sep = "."), width = width, height = height)
  plotIndex(df = slidwin, chr = chr, len = len, nSNPs = "nSNPs", n = minN, 
            index = "LB.index", band = 0.005, ylim = c(0, 1),
            ylab = paste(lowB, "SNP index", sep = " "), size = 1, axis.size = 1, 
            axis.lab.size = 1, axis.title.size = 1.5, type = "p")
  dev.off()
  png(filename = paste(outPrefix, lowB, "SNP_index.point.png", sep = "."), width = width, height = height, units = "in", res = 500)
  plotIndex(df = slidwin, chr = chr, len = len, nSNPs = "nSNPs", n = minN, 
            index = "LB.index", band = 0.005, ylim = c(0, 1),
            ylab = paste(lowB, "SNP index", sep = " "), size = 1, axis.size = 1, 
            axis.lab.size = 1, axis.title.size = 1.5, type = "p")
  dev.off()
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
  png(filename = paste(outPrefix, "ED4.point.png", sep = "."), width = width, height = height, units = "in", res = 500)
  plotIndex(df = slidwin, chr = chr, len = len, nSNPs = "nSNPs", n = minN, 
            index = "ED4", band = 0.005, ylim = c(0, max(slidwin$ED4, na.rm = TRUE)),
            ylab = "ED^4", size = 1, axis.size = 1, 
            axis.lab.size = 1, axis.title.size = 1.5, type = "p")
  dev.off()
  
  #### line plot
  ## Delta SNP index
  pdf(file = paste(outPrefix, "delta_SNP_index.line.pdf", sep = "."), width = width, height = height)
  plotIndex(df = slidwin, chr = chr, len = len, nSNPs = "nSNPs", n = minN, 
            index = "delta.index", band = 0.005, ylim = c(-1, 1),
            ylab = "Delta SNP index", size = 1.5, axis.size = 1, 
            axis.lab.size = 1, axis.title.size = 1.5, type = "l")
  dev.off()
  png(filename = paste(outPrefix, "delta_SNP_index.line.png", sep = "."), width = width, height = height, units = "in", res = 500)
  plotIndex(df = slidwin, chr = chr, len = len, nSNPs = "nSNPs", n = minN, 
            index = "delta.index", band = 0.005, ylim = c(-1, 1),
            ylab = "Delta SNP index", size = 1.5, axis.size = 1, 
            axis.lab.size = 1, axis.title.size = 1.5, type = "l")
  dev.off()
  ## High Bulk SNP index
  pdf(file = paste(outPrefix, highB, "SNP_index.line.pdf", sep = "."), width = width, height = height)
  plotIndex(df = slidwin, chr = chr, len = len, nSNPs = "nSNPs", n = minN, 
            index = "HB.index", band = 0.005, ylim = c(0, 1),
            ylab = paste(highB, "SNP index", sep = " "), size = 1.5, axis.size = 1, 
            axis.lab.size = 1, axis.title.size = 1.5, type = "l")
  dev.off()
  png(filename = paste(outPrefix, highB, "SNP_index.line.png", sep = "."), width = width, height = height, units = "in", res = 500)
  plotIndex(df = slidwin, chr = chr, len = len, nSNPs = "nSNPs", n = minN, 
            index = "HB.index", band = 0.005, ylim = c(0, 1),
            ylab = paste(highB, "SNP index", sep = " "), size = 1.5, axis.size = 1, 
            axis.lab.size = 1, axis.title.size = 1.5, type = "l")
  dev.off()
  ## Low Bulk SNP index
  pdf(file = paste(outPrefix, lowB, "SNP_index.line.pdf", sep = "."), width = width, height = height)
  plotIndex(df = slidwin, chr = chr, len = len, nSNPs = "nSNPs", n = minN, 
            index = "LB.index", band = 0.005, ylim = c(0, 1),
            ylab = paste(lowB, "SNP index", sep = " "), size = 1.5, axis.size = 1, 
            axis.lab.size = 1, axis.title.size = 1.5, type = "l")
  dev.off()
  png(filename = paste(outPrefix, lowB, "SNP_index.line.png", sep = "."), width = width, height = height, units = "in", res = 500)
  plotIndex(df = slidwin, chr = chr, len = len, nSNPs = "nSNPs", n = minN, 
            index = "LB.index", band = 0.005, ylim = c(0, 1),
            ylab = paste(lowB, "SNP index", sep = " "), size = 1.5, axis.size = 1, 
            axis.lab.size = 1, axis.title.size = 1.5, type = "l")
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
}else if (!is.na(CI)) {
  load(CI)
  dd6 <- dd6 %>% left_join(dltIndex_CI, by = c("HB.DP", "LB.DP"))
  slidwin <- winScan(x = dd6,
                     groups = "CHROM",
                     position = "POS",
                     values = c("HB.index", "LB.index", "delta.index", "ED", 
                                "CI95upper", "CI95lower", "CI99upper", "CI99lower"),
                     win_size = winSize,
                     win_step = winStep,
                     funs = c("mean")) %>% as_tibble() %>%
    dplyr::select(CHROM, win_start, win_end, win_mid, 
                  HB.index = HB.index_mean, LB.index = LB.index_mean, delta.index = delta.index_mean, 
                  CI95upper = CI95upper_mean, CI95lower = CI95lower_mean, 
                  CI99upper = CI99upper_mean, CI99lower = CI99lower_mean, 
                  ED= ED_mean, nSNPs = delta.index_n) %>%
    dplyr::mutate(ED4 = ED^4)
  write_tsv(x = slidwin, file = paste(outPrefix, "SlidingWindow.txt", sep = "."))
  write_csv(x = slidwin, file = paste(outPrefix, "SlidingWindow.csv", sep = "."))
  
  slidwin$POS = slidwin$win_mid
  #len <- chr %>% left_join(chr, by = "CHROM")
  #slidwin$CHROM <- factor(slidwin$CHROM, levels = chr$CHROM)
  ##
  pdf(file = paste(outPrefix, "delta_SNP_index.99CI.line.pdf", sep = "."), width = width, height = height)
  plotIndex(df = slidwin, chr = chr, len = len, CI = "CI99", nSNPs = "nSNPs", n = minN, 
            index = "delta.index", band = 0.005, #color = c("#177cb0","#f36838"), 
            ylim = c(-1, 1), ylab = "Delta SNP index", size = 1.5, axis.size = 1, 
            axis.lab.size = 1, axis.title.size = 1.5, type = "l")
  dev.off()
  png(filename = paste(outPrefix, "delta_SNP_index.99CI.line.png", sep = "."), width = width, height = height, units = "in", res = 500)
  plotIndex(df = slidwin, chr = chr, len = len, CI = "CI99", nSNPs = "nSNPs", n = minN, 
            index = "delta.index", band = 0.005, #color = c("#177cb0","#f36838"), 
            ylim = c(-1, 1), ylab = "Delta SNP index", size = 1.5, axis.size = 1, 
            axis.lab.size = 1, axis.title.size = 1.5, type = "l")
  dev.off()
  ##
  pdf(file = paste(outPrefix, "delta_SNP_index.95CI.line.pdf", sep = "."), width = width, height = height)
  plotIndex(df = slidwin, chr = chr, len = len, CI = "CI95", nSNPs = "nSNPs", n = minN, 
            index = "delta.index", band = 0.005, #color = c("#177cb0","#f36838"), 
            ylim = c(-1, 1), ylab = "Delta SNP index", size = 1.5, axis.size = 1, 
            axis.lab.size = 1, axis.title.size = 1.5, type = "l")
  dev.off()
  png(filename = paste(outPrefix, "delta_SNP_index.95CI.line.png", sep = "."), width = width, height = height, units = "in", res = 500)
  plotIndex(df = slidwin, chr = chr, len = len, CI = "CI95", nSNPs = "nSNPs", n = minN, 
            index = "delta.index", band = 0.005, #color = c("#177cb0","#f36838"), 
            ylim = c(-1, 1), ylab = "Delta SNP index", size = 1.5, axis.size = 1, 
            axis.lab.size = 1, axis.title.size = 1.5, type = "l")
  dev.off()
  ##
  pdf(file = paste(outPrefix, highB, "SNP_index.line.pdf", sep = "."), width = width, height = height)
  plotIndex(df = slidwin, chr = chr, len = len, nSNPs = "nSNPs", n = minN, 
            index = "HB.index", band = 0.005, #color = c("#177cb0","#f36838"), 
            ylim = c(0, 1), ylab = paste(highB, "SNP index", sep = " "), size = 1.5, axis.size = 1, 
            axis.lab.size = 1, axis.title.size = 1.5, type = "l")
  dev.off()
  png(filename = paste(outPrefix, highB, "SNP_index.line.png", sep = "."), width = width, height = height, units = "in", res = 500)
  plotIndex(df = slidwin, chr = chr, len = len, nSNPs = "nSNPs", n = minN, 
            index = "HB.index", band = 0.005, #color = c("#177cb0","#f36838"), 
            ylim = c(0, 1), ylab = paste(highB, "SNP index", sep = " "), size = 1.5, axis.size = 1, 
            axis.lab.size = 1, axis.title.size = 1.5, type = "l")
  dev.off()
  ##
  pdf(file = paste(outPrefix, lowB, "SNP_index.line.pdf", sep = "."), width = width, height = height)
  plotIndex(df = slidwin, chr = chr, len = len, nSNPs = "nSNPs", n = minN, 
            index = "LB.index", band = 0.005, #color = c("#177cb0","#f36838"), 
            ylim = c(0, 1), ylab = paste(lowB, "SNP index", sep = " "), size = 1.5, axis.size = 1, 
            axis.lab.size = 1, axis.title.size = 1.5, type = "l")
  dev.off()
  png(filename = paste(outPrefix, highB, "SNP_index.line.png", sep = "."), width = width, height = height, units = "in", res = 500)
  plotIndex(df = slidwin, chr = chr, len = len, nSNPs = "nSNPs", n = minN, 
            index = "LB.index", band = 0.005, #color = c("#177cb0","#f36838"), 
            ylim = c(0, 1), ylab = paste(lowB, "SNP index", sep = " "), size = 1.5, axis.size = 1, 
            axis.lab.size = 1, axis.title.size = 1.5, type = "l")
  dev.off()
  ##
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
  ##
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
  
  #### point plot
  ##
  pdf(file = paste(outPrefix, "delta_SNP_index.99CI.point.pdf", sep = "."), width = width, height = height)
  plotIndex(df = slidwin, chr = chr, len = len, CI = "CI99", nSNPs = "nSNPs", n = minN, 
            index = "delta.index", band = 0.005, #color = c("#177cb0","#f36838"), 
            ylim = c(-1, 1), ylab = "Delta SNP index", size = 1, axis.size = 1, 
            axis.lab.size = 1, axis.title.size = 1.5, type = "p")
  dev.off()
  png(filename = paste(outPrefix, "delta_SNP_index.99CI.point.png", sep = "."), width = width, height = height, units = "in", res = 500)
  plotIndex(df = slidwin, chr = chr, len = len, CI = "CI99", nSNPs = "nSNPs", n = minN, 
            index = "delta.index", band = 0.005, #color = c("#177cb0","#f36838"), 
            ylim = c(-1, 1), ylab = "Delta SNP index", size = 1, axis.size = 1, 
            axis.lab.size = 1, axis.title.size = 1.5, type = "p")
  dev.off()
  ##
  pdf(file = paste(outPrefix, "delta_SNP_index.95CI.point.pdf", sep = "."), width = width, height = height)
  plotIndex(df = slidwin, chr = chr, len = len, CI = "CI95", nSNPs = "nSNPs", n = minN, 
            index = "delta.index", band = 0.005, #color = c("#177cb0","#f36838"), 
            ylim = c(-1, 1), ylab = "Delta SNP index", size = 1, axis.size = 1, 
            axis.lab.size = 1, axis.title.size = 1.5, type = "p")
  dev.off()
  png(filename = paste(outPrefix, "delta_SNP_index.95CI.point.png", sep = "."), width = width, height = height, units = "in", res = 500)
  plotIndex(df = slidwin, chr = chr, len = len, CI = "CI95", nSNPs = "nSNPs", n = minN, 
            index = "delta.index", band = 0.005, #color = c("#177cb0","#f36838"), 
            ylim = c(-1, 1), ylab = "Delta SNP index", size = 1, axis.size = 1, 
            axis.lab.size = 1, axis.title.size = 1.5, type = "p")
  dev.off()
  ##
  pdf(file = paste(outPrefix, highB, "SNP_index.point.pdf", sep = "."), width = width, height = height)
  plotIndex(df = slidwin, chr = chr, len = len, nSNPs = "nSNPs", n = minN, 
            index = "HB.index", band = 0.005, #color = c("#177cb0","#f36838"), 
            ylim = c(0, 1), ylab = paste(highB, "SNP index", sep = " "), size = 1, axis.size = 1, 
            axis.lab.size = 1, axis.title.size = 1.5, type = "p")
  dev.off()
  png(filename = paste(outPrefix, highB, "SNP_index.point.png", sep = "."), width = width, height = height, units = "in", res = 500)
  plotIndex(df = slidwin, chr = chr, len = len, nSNPs = "nSNPs", n = minN, 
            index = "HB.index", band = 0.005, #color = c("#177cb0","#f36838"), 
            ylim = c(0, 1), ylab = paste(highB, "SNP index", sep = " "), size = 1, axis.size = 1, 
            axis.lab.size = 1, axis.title.size = 1.5, type = "p")
  dev.off()
  ##
  pdf(file = paste(outPrefix, lowB, "SNP_index.point.pdf", sep = "."), width = width, height = height)
  plotIndex(df = slidwin, chr = chr, len = len, nSNPs = "nSNPs", n = minN, 
            index = "LB.index", band = 0.005, #color = c("#177cb0","#f36838"), 
            ylim = c(0, 1), ylab = paste(lowB, "SNP index", sep = " "), size = 1, axis.size = 1, 
            axis.lab.size = 1, axis.title.size = 1.5, type = "p")
  dev.off()
  png(filename = paste(outPrefix, lowB, "SNP_index.point.png", sep = "."), width = width, height = height, units = "in", res = 500)
  plotIndex(df = slidwin, chr = chr, len = len, nSNPs = "nSNPs", n = minN, 
            index = "LB.index", band = 0.005, #color = c("#177cb0","#f36838"), 
            ylim = c(0, 1), ylab = paste(lowB, "SNP index", sep = " "), size = 1, axis.size = 1, 
            axis.lab.size = 1, axis.title.size = 1.5, type = "p")
  dev.off()
  ##
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
  ##
  pdf(file = paste(outPrefix, "ED4.point.pdf", sep = "."), width = width, height = height)
  plotIndex(df = slidwin, chr = chr, len = len, nSNPs = "nSNPs", n = minN, 
            index = "ED4", band = 0.005, ylim = c(0, max(slidwin$ED4, na.rm = TRUE)),
            ylab = "ED^4", size = 1, axis.size = 1, 
            axis.lab.size = 1, axis.title.size = 1.5, type = "p")
  dev.off()
  png(filename = paste(outPrefix, "ED4.point.png", sep = "."), width = width, height = height, units = "in", res = 500)
  plotIndex(df = slidwin, chr = chr, len = len, nSNPs = "nSNPs", n = minN, 
            index = "ED4", band = 0.005, ylim = c(0, max(slidwin$ED4, na.rm = TRUE)),
            ylab = "ED^4", size = 1, axis.size = 1, 
            axis.lab.size = 1, axis.title.size = 1.5, type = "p")
  dev.off()
  
}


CI95 <- getQTL(data = slidwin, CI = 95, n = minN, export = TRUE, filename = paste(outPrefix, "95CI.csv", sep = "."))
CI99 <- getQTL(data = slidwin, CI = 99, n = minN, export = TRUE, filename = paste(outPrefix, "99CI.csv", sep = "."))


## plot target chrom
options(scipen=200)
plotTargetChrom(df = slidwin, CI = 95, minN = minN, outPrefix = outPrefix, chr = chr, len = len)
plotTargetChrom(df = slidwin, CI = 99, minN = minN, outPrefix = outPrefix, chr = chr, len = len)
