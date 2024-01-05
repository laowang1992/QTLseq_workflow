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

sample <- read_tsv(file = sampleInfo, col_names = F, show_col_types = FALSE) %>% pull(X1)
chr <- read_tsv(file = chrInfo, col_names = c("CHROM", "LABEL"), show_col_types = FALSE)

chromInfo <- read_tsv(file = chrLen, col_names = F, show_col_types = FALSE) %>% 
  select(CHROM = X1, LEN = X2) %>% 
  right_join(chr, by = "CHROM") %>% 
  mutate(start = 0) %>% 
  select(chr = LABEL, start, end = LEN)

mid <- 0
for (i in seq_along(sample)) {
  depth <- read_tsv(file = paste(sample[i], "dd.window.depth", sep = "."), col_names = F, show_col_types = FALSE)
  if (mid < median(depth$X4)) {
    mid <- median(depth$X4)
  }
}

getPalette = colorRampPalette(brewer.pal(8, "Set2"))
color = getPalette(length(sample))

n <- length(sample)

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
  for(i in seq_along(sample)){
    cat(date(), ", read and plot ", sample[i], " ...\n", sep = "")
    depth <- read_tsv(file = paste(sample[i], "dd.window.depth", sep = "."), col_names = F, show_col_types = FALSE) %>%
      rename(CHROM = X1, start = X2, end = X3, value = X4) %>%
      right_join(chr, by = "CHROM") %>% select(chr = LABEL, start, end, value) %>%
      mutate(pos = (start+end)/2)
    p <- ggplot(depth, aes(x = pos, y = value, fill = chr)) +
      geom_area() +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, 2.5*mid), n.breaks = 2) +
      scale_fill_manual(values = getPalette(nrow(chromInfo))) +
      labs(x = NULL, y = "Depth") +
      facet_grid(chr ~ .) +
      cowplot::theme_half_open() +
      theme(strip.text.y = element_text(angle = 0),
            strip.background = element_rect(color = NA, fill = NA),
            legend.position = "NULL")
    ggsave(p, filename = paste(sample[i], "CoverageDepth.pdf", sep = "."), width = 8, height = dim(chr)[[1]] * 0.45 + 0.5)
    ggsave(p, filename = paste(sample[i], "CoverageDepth.png", sep = "."), width = 8, height = dim(chr)[[1]] * 0.45 + 0.5, dpi = 500)
  }
  stop("Too many samples, draw depth respectively.")
}

png(filename = "CoverageDepth.png", width = width, height = height, units = "in", res = 500)
par(mar = mar + 0.1)
#par(mar = c(5, 4, 4, 2) + 0.1)
circos.par("start.degree" = 90, track.height = track.height, track.margin = c(0, 0), cell.padding = c(0, 1.00, 0.02, 1.00))
cat(date(), ", initialization ...\n", sep = "")
circos.genomicInitialize(chromInfo, plotType = c("axis", "labels"))

for (i in seq_along(sample)) {
  cat(date(), ", read and plot ", sample[i], " ...\n", sep = "")
  depth <- read_tsv(file = paste(sample[i], "dd.window.depth", sep = "."), col_names = F, show_col_types = FALSE) %>%
    rename(CHROM = X1, start = X2, end = X3, value = X4) %>%
    right_join(chr, by = "CHROM") %>% select(chr = LABEL, start, end, value) %>%
    mutate(value = if_else(value > 2.5*mid, 2.5*mid, value))
  circos.genomicTrack(depth, ylim = c(0, 2.5*mid), bg.border = NA,
                      panel.fun = function(region, value, ...) {
                        circos.genomicLines(region, value, area = TRUE, col = color[i], border = NA)
                      })
}

circos.clear()
legend(position, inset=.05, sample, pch=15, col=color, pt.cex = pt.cex, cex = cex, box.col = NA)
dev.off()

# 
pdf(file = "CoverageDepth.pdf", width = width, height = height)
par(mar = mar + 0.1)
#par(mar = c(5, 4, 4, 2) + 0.1)
circos.par("start.degree" = 90, track.height = track.height, track.margin = c(0, 0), cell.padding = c(0, 1.00, 0.02, 1.00))
cat(date(), ", initialization ...\n", sep = "")
circos.genomicInitialize(chromInfo, plotType = c("axis", "labels"))

for (i in seq_along(sample)) {
  cat(date(), ", read and plot ", sample[i], " ...\n", sep = "")
  depth <- read_tsv(file = paste(sample[i], "dd.window.depth", sep = "."), col_names = F, show_col_types = FALSE) %>%
    rename(CHROM = X1, start = X2, end = X3, value = X4) %>%
    right_join(chr, by = "CHROM") %>% select(chr = LABEL, start, end, value) %>%
    mutate(value = if_else(value > 2.5*mid, 2.5*mid, value))
  circos.genomicTrack(depth, ylim = c(0, 2.5*mid), bg.border = NA,
                      panel.fun = function(region, value, ...) {
                        circos.genomicLines(region, value, area = TRUE, col = color[i], border = NA)
                      })
}

circos.clear()
legend(position, inset=.05, sample, pch=15, col=color, pt.cex = pt.cex, cex = cex, box.col = NA)
dev.off()

