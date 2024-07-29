#!/usr/bin/env Rscript
# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
# Create a parser
p <- arg_parser("Coverage Depth Circos plot")

## Add command line arguments
#
p <- add_argument(p, "--sampleInfo", help = "Sample list file, tab seperated with first colume showing samples", type = "character")
p <- add_argument(p, "--chrInfo", help = "A file containing chr included in circos and LABEL, two columes (CHROM\tLABEL), tab seperated", type = "character")
p <- add_argument(p, "--chrLen", help = "A file containing chr length", type = "character")
# Parse the command line arguments
argv <- parse_args(p)

sampleInfo <- argv$sampleInfo
chrInfo <- argv$chrInfo
chrLen <- argv$chrLen

library(circlize)
library(tidyverse)
library(RColorBrewer)

test <- FALSE
if (test) {
  sampleInfo <- "./samples.txt"
  chrInfo <- "./chrom.txt"
  chrLen <- "./ref.len"
}


samples <- read_tsv(file = sampleInfo, col_names = F, show_col_types = FALSE) %>% pull(X1)
# 检查一下有没有重复样本或没有样本
if (length(samples) != length(unique(samples)) || length(samples) == 0) {
  stop("Duplication samples or null sample.\n")
}

chr <- read_tsv(file = chrInfo, col_names = c("Chr", "Name"), col_types = cols(Chr = "c"), show_col_types = FALSE)
chromInfo <- read_tsv(file = chrLen, col_names = c("Chr", "End"), col_types = cols(Chr = "c"), show_col_types = FALSE) %>% 
  right_join(chr, by = "Chr") %>% 
  mutate(Start = 0) %>% 
  select(Chr, Name, Start, End)

# median
mid <- 0
for (sample in samples) {
  subdf <- read_tsv(file = paste(sample, "win.stat.gz", sep = "."), col_names = T, col_types = cols(`#Chr` = "c"), show_col_types = FALSE, comment = "##") %>% 
    mutate(Sample = sample) %>% rename(Chr = `#Chr`)
  subdf <- chr %>% select(Chr) %>% left_join(subdf, by = "Chr")
  if (mid < median(subdf$MeanDepth)) {
    mid <- median(subdf$MeanDepth)
  }
  if (match(sample, samples) == 1) {
    df <- subdf
  } else if (match(sample, samples) > 1) {
    df <- rbind(df, subdf)
  }
}

n <- length(samples)

getPalette = colorRampPalette(brewer.pal(8, "Set2"))
color = getPalette(n)

# c(bottom, left, top, right)
if (n <= 6 & n > 0) {
  width = 6
  height = 6
  mar = c(1, 1, 1, 1)
  track.height = 0.1
  position = "center"
  pt.cex = 1.5
  cex = 0.8
} else if (n > 6 & n <=10) {
  width = 7
  height = 6
  mar = c(1, 1, 1, 1)
  track.height = 0.7 / n
  position = "topright"
  pt.cex = 1.1
  cex = 0.6
} else if (n > 10){
  for(sample in samples){
    i <- match(sample, samples)
    cat(date(), ", read and plot ", sample, " ...\n", sep = "")
    subdf <- df %>% filter(Sample == sample) %>% right_join(chr, by = "Chr") %>% 
      mutate(Pos = as.integer((Start+End)/2))
    p <- ggplot(subdf, aes(x = Pos, y = MeanDepth, fill = Chr)) +
      geom_area() +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(n.breaks = 2) +
      scale_fill_manual(values = getPalette(nrow(chromInfo))) +
      labs(x = NULL, y = "Depth") + 
      coord_cartesian(ylim = c(0, 2.5*mid)) + 
      facet_grid(Name ~ .) + 
      cowplot::theme_half_open() + 
      theme(strip.text.y = element_text(angle = 0),
            strip.background = element_rect(color = NA, fill = NA),
            legend.position = "NULL")
    ggsave(p, filename = paste(sample, "CoverageDepth.pdf", sep = "."), width = 8, height = dim(chr)[[1]] * 0.45 + 0.5)
    ggsave(p, filename = paste(sample, "CoverageDepth.png", sep = "."), width = 8, height = dim(chr)[[1]] * 0.45 + 0.5, dpi = 500)
  }
  stop("Too many samples, draw depth respectively.")
}

png(filename = "CoverageDepth.png", width = width, height = height, units = "in", res = 500)
par(mar = mar + 0.1)
#par(mar = c(5, 4, 4, 2) + 0.1)
circos.par("start.degree" = 90, track.height = track.height, track.margin = c(0, 0), cell.padding = c(0, 1.00, 0.02, 1.00))
cat(date(), ", initialization ...\n", sep = "")
circos.genomicInitialize(chromInfo %>% select(-Chr), plotType = c("axis", "labels"))

for (sample in samples) {
  i <- match(sample, samples)
  cat(date(), ", read and plot ", sample, " ...\n", sep = "")
  subdf <- df %>% filter(Sample == sample) %>% right_join(chr, by = "Chr") %>% 
    select(Name, Start, End, MeanDepth) %>% 
    mutate(MeanDepth = if_else(MeanDepth > 2.5*mid, 2.5*mid, MeanDepth))
  circos.genomicTrack(subdf, ylim = c(0, 2.5*mid), bg.border = NA,
                      panel.fun = function(region, value, ...) {
                        circos.genomicLines(region, value, area = TRUE, col = color[i], border = NA)
                      })
}
circos.clear()
legend(position, inset=.05, samples, pch=15, col=color, pt.cex = pt.cex, cex = cex, box.col = NA)
dev.off()

pdf(file = "CoverageDepth.pdf", width = width, height = height)
par(mar = mar + 0.1)
#par(mar = c(5, 4, 4, 2) + 0.1)
circos.par("start.degree" = 90, track.height = track.height, track.margin = c(0, 0), cell.padding = c(0, 1.00, 0.02, 1.00))
cat(date(), ", initialization ...\n", sep = "")
circos.genomicInitialize(chromInfo %>% select(-Chr), plotType = c("axis", "labels"))

for (sample in samples) {
  i <- match(sample, samples)
  cat(date(), ", read and plot ", sample, " ...\n", sep = "")
  subdf <- df %>% filter(Sample == sample) %>% right_join(chr, by = "Chr") %>% 
    select(Name, Start, End, MeanDepth) %>% 
    mutate(MeanDepth = if_else(MeanDepth > 2.5*mid, 2.5*mid, MeanDepth))
  circos.genomicTrack(subdf, ylim = c(0, 2.5*mid), bg.border = NA,
                      panel.fun = function(region, value, ...) {
                        circos.genomicLines(region, value, area = TRUE, col = color[i], border = NA)
                      })
}
circos.clear()
legend(position, inset=.05, samples, pch=15, col=color, pt.cex = pt.cex, cex = cex, box.col = NA)
dev.off()

# coverage rate
covRate <- df %>% group_by(Sample) %>% summarise(TotalLen=sum(Length), CoveredSite=sum(CoveredSite)) %>% 
  mutate(CovRate = CoveredSite/TotalLen) %>% select(-TotalLen)
write_tsv(x = covRate, file = "cov_stat.txt")
write_csv(x = covRate, file = "cov_stat.csv")
