#!/usr/bin/env Rscript
# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
# Create a parser
p <- arg_parser("a program for QTLseq, TABEL file from GATK as input")

## Add command line arguments
#
p <- add_argument(p, "--list", short = "-l", help = "list table header", flag = TRUE)
p <- add_argument(p, "--createParameter", short = "-C", help = "create a Parameter csv file", flag = TRUE)

p <- add_argument(p, "--parameter", help = "Parameter csv file containing parameters, if the same parameter is listed in CMD line, this CMD line parameter will be omitted", type = "character")
#
p <- add_argument(p, "--input", help = "TABEL file from GATK", type = "character")
p <- add_argument(p, "--out", help = "A prefix for output file", type = "character")
#
p <- add_argument(p, "--CI", help = "QTLseqCI**.RData, you must check whether it meet with you sample, e.g.: popType, bulk size, depth", type = "character", default = NULL)
#
p <- add_argument(p, "--highP", help = "The parent name with high phenotype", type = "character")
p <- add_argument(p, "--lowP", help = "The parent name with low phenotype", type = "character")
p <- add_argument(p, "--highB", help = "The bulk name with high phenotype", type = "character")
p <- add_argument(p, "--lowB", help = "The bulk name with low phenotype", type = "character")
p <- add_argument(p, "--popType", help = "Population type, 'F2' or 'RIL'", type = "character")
# minQ这个参数其实没有用，但是如果删除，运行以前的命令行会报错，所以现在先这样
p <- add_argument(p, "--minQ", help = "Minimum QUAL for variation", type = "numeric", default = 30)
#
p <- add_argument(p, "--bulkSizeH", help = "Bulk size with high phenotype", type = "numeric")
p <- add_argument(p, "--bulkSizeL", help = "Bulk size with low phenotype", type = "numeric")

p <- add_argument(p, "--minGQ", help = "Minimum GQ for parent genotype", type = "numeric", default = 20)
#
p <- add_argument(p, "--minHPdp", help = "Minimum depth for high parent", type = "numeric", default = 6)
p <- add_argument(p, "--minLPdp", help = "Minimum depth for low parent", type = "numeric", default = 6)
p <- add_argument(p, "--minHBdp", help = "Minimum depth for high bulk", type = "numeric", default = 6)
p <- add_argument(p, "--minLBdp", help = "Minimum depth for low bulk", type = "numeric", default = 6)
# 下面四个参数其实也没有用，真正的最大允许深度根据ave+3*sd计算得出，只是为了兼容以前的命令行参数
p <- add_argument(p, "--maxHPdp", help = "Maxmum depth for high parent", type = "numeric", default = 100)
p <- add_argument(p, "--maxLPdp", help = "Maxmum depth for low parent", type = "numeric", default = 100)
p <- add_argument(p, "--maxHBdp", help = "Maxmum depth for high bulk", type = "numeric", default = 100)
p <- add_argument(p, "--maxLBdp", help = "Maxmum depth for low bulk", type = "numeric", default = 100)

#
p <- add_argument(p, "--winSize", help = "Window size for sliding window statistics", type = "numeric", default = 2000000)
p <- add_argument(p, "--winStep", help = "Window step for sliding window statistics", type = "numeric", default = 200000)
#
p <- add_argument(p, "--minN", help = "Minimum SNP number in a window for plot", type = "numeric", default = 20)
#
p <- add_argument(p, "--width", help = "Delta SNP index plot width", type = "numeric", default = 12)
p <- add_argument(p, "--height", help = "Delta SNP index plot height", type = "numeric", default = 4)

# Parse the command line arguments
argv <- parse_args(p)
filename <- argv$input

#
if (argv$list) {
  cat("This is the header row:\n")
  df <- read.table(argv$input, nrows = 1, header = T)
  col_name <- colnames(df)
  print(col_name)
  quit(save = "no")
}
# 创建一个参数模板文件
if (argv$createParameter) {
  cat("You will get a template file of parameter.\n")
  parameter <- data.frame(outPrefix = argv$out, CI = argv$CI, 
                          highParent = argv$highP, lowParent = argv$lowP, highBulk = argv$highB, lowBulk = argv$lowB, 
                          populationType = argv$popType, highBulkSize = argv$bulkSizeH, lowBulkSize = argv$bulkSizeL, 
                          minGQ = argv$minGQ, 
                          minHPdp = argv$minHPdp, minLPdp = argv$minLPdp, 
                          minHBdp = argv$minHBdp, minLBdp = argv$minLBdp, 
                          winSize = argv$winSize, winStep = argv$winStep, minN = argv$minN, 
                          width = argv$width, height = argv$height
                          )
  write.csv(x = parameter, file = "Parameter.csv", na = "", row.names = F)
  quit(save = "no")
}

library(tidyverse)
#library(QTLseqr)
#library(windowscanr)
library(cowplot)
#library(CMplot)
library(ggsci)
library(RColorBrewer)


if (FALSE) {
  ## file name
  filename <- "./BSA.filter.SNPs.table.gz"
  outPrefix <- "RY"
  
  CI <- "./QTLseqCI.indH30_indL30_F2.Depth_5_200.Rep_10000.RData"
  ##
  bulkSizeH <- 30
  bulkSizeL <- 30
  ## sample name R3	qY	R	Y
  highP <- "R3"
  lowP <- "qY"
  highB <- "R"
  lowB <- "Y"
  ##
  popType <- "F2"  # F2 ro RIL
  ##
  minGQ <- 20
  ## filter depth parameter, only set min depth, max is ave+3*sd
  minHPdp <- 6
  minLPdp <- 6
  minHBdp <- 6
  minLBdp <- 6

  ## sliding window parameter
  winSize <- 2000000
  winStep <- 200000
  #
  minN <- 50
  ##
  width <- 12
  height <- 4
}

# 读取数据和函数
chr <- read_tsv("./chrom.txt", show_col_types = FALSE)
len <- read_tsv(file = "./ref.len", col_names = c("CHROM", "Len"), show_col_types = FALSE)
source("./Support_functions.R")
# 读取基因型等相关信息
data <- read_tsv(file = filename, show_col_types = FALSE)
# 当前目录
wd <- getwd()
# 
if (is.na(argv$parameter)) {
  cat("Get parameters from CMD line.\n")
  parameter <- data.frame(filename = argv$input, outPrefix = argv$out, CI = argv$CI, minGQ = argv$minGQ, 
                          highParent = argv$highP, lowParent = argv$lowP, highBulk = argv$highB, lowBulk = argv$lowB, 
                          populationType = argv$popType, highBulkSize = argv$bulkSizeH, lowBulkSize = argv$bulkSizeL, 
                          minHPdp = argv$minHPdp, minLPdp = argv$minLPdp, minHBdp = argv$minHBdp, minLBdp = argv$minLBdp, 
                          winSize = argv$winSize, winStep = argv$winStep, minN = argv$minN, 
                          width = argv$width, height = argv$height
                          )
} else {
  cat("Get parameters from parameter csv file.\n")
  parameter <- read.csv(argv$parameter, header = T)
}

for (i in seq_along(parameter[,1])) {
  outPrefix <- parameter$outPrefix[i]
  
  CI <- parameter$CI[i]
  #
  bulkSizeH <- parameter$highBulkSize[i]
  bulkSizeL <- parameter$lowBulkSize[i]
  # sample name
  highP <- parameter$highParent[i]
  lowP <- parameter$lowParent[i]
  highB <- parameter$highBulk[i]
  lowB <- parameter$lowBulk[i]
  #
  popType <- parameter$populationType[i]  # F2 ro RIL
  #
  minGQ <- parameter$minGQ[i]
  # filter depth parameter
  minHPdp <- parameter$minHPdp[i]
  minLPdp <- parameter$minLPdp[i]
  minHBdp <- parameter$minHBdp[i]
  minLBdp <- parameter$minLBdp[i]
  
  # sliding window parameter
  winSize <- parameter$winSize[i]
  winStep <- parameter$winStep[i]
  #
  minN <- parameter$minN[i]
  #
  width <- parameter$width[i]
  height <- parameter$height[i]
  
  if (!is.na(highP) && paste(highP, "GT", sep = ".") %in% colnames(data) && 
      !is.na(lowP) && paste(lowP, "GT", sep = ".") %in% colnames(data)) {
    cat("There if 2 parents\n")
    source("./for_2_parent.R")
  } else if (!is.na(highP) && paste(highP, "GT", sep = ".") %in% colnames(data)) {
    cat("There is only", highP, "parent\n")
    source("./for_1_parent.R")
  } else if (!is.na(lowP) && paste(lowP, "GT", sep = ".") %in% colnames(data)) {
    cat("There is only", lowP, "parent\n")
    tmp <- highP; highP <- lowP; lowP <- tmp
    tmp <- highB; highB <- lowB; lowB <- tmp
    tmp <- bulkSizeH; bulkSizeH <- bulkSizeL; bulkSizeL <- tmp
    tmp <- minHPdp; minHPdp <- minLPdp; minLPdp <- tmp
    tmp <- minHBdp; minHBdp <- minLBdp; minLBdp <- tmp
    source("./for_1_parent.R")
  } else {
    cat("There is 0 parent\n")
    source("./for_0_parent.R")
  }
  # 目录设置
  dir.create(path = outPrefix, showWarnings = FALSE)
  # 这里写的不是很顺滑，以后应该考虑修改
  if (!is.na(CI)) {
    file.copy(from = CI, to = paste(outPrefix, CI, sep = "/"))
  }
  setwd(outPrefix)
  write.table(x = parameter[i,], paste(outPrefix, "parameter.txt", sep = "."), quote = F, sep = "\t", row.names = F)
  df <- selectSample_and_filterGT(data = data, highP = highP, lowP = lowP, highB = highB, lowB = lowB, minGQ = minGQ)
  depth_statistics(data = df, outPrefix = outPrefix, highP = highP, lowP = lowP, highB = highB, lowB = lowB)
  final_dp <- filterDP(data = df, minHPdp = minHPdp, minLPdp = minLPdp, minHBdp = minHBdp, minLBdp = minLBdp)
  SNP_distribution(data = final_dp, outPrefix = outPrefix, chr = chr)
  export_dp(data = final_dp, outPrefix = outPrefix, highP = highP, lowP = lowP, highB = highB, lowB = lowB)
  slidwin <- calc_index_etc(data = final_dp, CI = CI, outPrefix = outPrefix, 
                            minHBdp = minHBdp, minLBdp = minLBdp, popType = popType, 
                            bulkSizeH = bulkSizeH, bulkSizeL = bulkSizeL, 
                            winSize = winSize, winStep = winStep)
  export_figure(data = slidwin, outPrefix = outPrefix, chr = chr, len = len, minN = minN, highB = highB, lowB = lowB, width = width, height = height)
  # 如果提供至少一个亲本就可以计算区间
  if ((!is.na(highP) && paste(highP, "GT", sep = ".") %in% colnames(data)) || (!is.na(lowP) && paste(lowP, "GT", sep = ".") %in% colnames(data))) {
    getQTL_and_exportFigure(data = slidwin, outPrefix = outPrefix, minN = minN, chr = chr, len = len)
  }
  # 回到原始工作目录
  setwd(wd)
}

