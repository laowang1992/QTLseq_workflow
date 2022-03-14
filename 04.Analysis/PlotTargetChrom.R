library(tidyverse)
library(cowplot)
options(scipen=200)

inputPrefix <- "BSA"
CI <- "99"
qtlfile <- paste(inputPrefix, paste(CI, "CI", sep = ""), "csv", sep = ".")
indexfile <- paste(inputPrefix, "QTLseqrdeltaSNPindex", "csv", sep = ".")
minN <- 0

#
snp_index <- 
qtl <- read_csv(qtlfile) %>% 
  dplyr::filter(length > 0, nSNPs > minN)
index <- read_csv(indexfile)

for (i in seq_along(qtl$CHROM)) {
  chr <- qtl[i, "CHROM"][[1]]
  start <- qtl[i, "start"][[1]]
  end <- qtl[i, "end"][[1]]
  df <- index %>% filter(CHROM == chr)
  p <- ggplot(df) +
    geom_rect(aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),fill = '#FF3300', color = "#FF3300", alpha = .005) +
    geom_line(aes(x = POS, y = tricubeDeltaSNP), color = "blue") +
    geom_line(aes(x = POS, y = CI_99), color = "gray60") + 
    geom_line(aes(x = POS, y = -CI_99), color = "gray60") +
    geom_hline(aes(yintercept=0)) + 
    labs(x = chr, y = "delta SNP index") +
    coord_cartesian(ylim = c(-1, 1)) +
    theme_half_open()
  ggsave(p, filename = paste(paste(chr, start, end, sep = "_"), paste(CI, "CI", sep = ""), "png", sep = "."), height = 3.5, width = 4.5, dpi = 500)
  print(paste(chr, "has been done...", sep = " "))
}

