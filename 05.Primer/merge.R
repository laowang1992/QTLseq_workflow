library(tidyverse)
#library(xlsx)

argv <- commandArgs(trailingOnly = T)

test <- FALSE
if (test) {
  argv[1] <- "public.indels.txt"
  argv[2] <- "public.primer.txt"
}

gt <- read_tsv(file = argv[1]) %>%
  mutate(REF_len = str_length(REF), ALT_len = str_length(ALT), 
         DIFF = abs(REF_len - ALT_len)) %>% distinct()
primer <- read_tsv(file = argv[2]) %>% distinct()

p0 <- primer %>% filter(NUM == 0) %>% select(ID, PrimerL1 = PrimerL, PrimerR1 = PrimerR, Hit1 = Hit, TmL1 = TmL, TmR1 = TmR, GCL1 = GCL, GCR1 = GCR, Length1 = Length)
p1 <- primer %>% filter(NUM == 1) %>% select(ID, PrimerL2 = PrimerL, PrimerR2 = PrimerR, Hit2 = Hit, TmL2 = TmL, TmR2 = TmR, GCL2 = GCL, GCR2 = GCR, Length2 = Length)
p2 <- primer %>% filter(NUM == 2) %>% select(ID, PrimerL3 = PrimerL, PrimerR3 = PrimerR, Hit3 = Hit, TmL3 = TmL, TmR3 = TmR, GCL3 = GCL, GCR3 = GCR, Length3 = Length)

df <- gt %>% left_join(p0, by = "ID") %>% 
  left_join(p1, by = "ID") %>%
  left_join(p2, by = "ID")


write_csv(x = df, file = str_replace(argv[2], "primer.txt", "primer.csv"))
#write.xlsx(x = df, file = str_replace(argv[2], "primer.txt", "primer.xlsx"))

# Primer distribution
library(circlize)
library(RColorBrewer)
source("../slidingWindow.R")

chrInfo <- "../../refseq/chrom.txt"
chrLen <- "../../refseq//ref.len"
chr <- read_tsv(file = chrInfo, col_names = c("CHROM", "Name"), col_types = cols(Chr = "c"), show_col_types = FALSE)
chromInfo <- read_tsv(file = chrLen, col_names = c("CHROM", "End"), col_types = cols(CHROM = "c"), show_col_types = FALSE) %>% 
  right_join(chr, by = "CHROM") %>% 
  mutate(Start = 0) %>% 
  select(CHROM, Name, Start, End)

df <- df %>% 
  mutate(havePrimer = if_else(is.na(PrimerL1) & is.na(PrimerL2) & is.na(PrimerL3), 0, 1), 
         isUnique   = if_else(Hit1 == 1 | Hit2 == 1 | Hit3 == 1, 1, 0))
primerDensity <- slidingWindow(df = df, winSize = 200000, winStep = 200000, groups = "CHROM", position = "POS", values = c("havePrimer", "isUnique"), fun = "sum")

width = 6
height = 6
mar = c(1, 1, 1, 1)
track.height = 0.1
position = "center"
pt.cex = 1.5
cex = 0.8

png(filename = str_replace(argv[2], "primer.txt", "Primer_Distribution.png"), width = width, height = height, units = "in", res = 500)
par(mar = mar + 0.1)
#par(mar = c(5, 4, 4, 2) + 0.1)
circos.par("start.degree" = 90, track.height = track.height, track.margin = c(0, 0), cell.padding = c(0, 1.00, 0.02, 1.00))
cat(date(), ", initialization ...\n", sep = "")
circos.genomicInitialize(chromInfo %>% select(-CHROM), plotType = c("axis", "labels"))

havePrimer <- chr %>% left_join(primerDensity, by = "CHROM") %>% select(Name, win_start, win_end, havePrimer_sum)
isUnique <- chr %>% left_join(primerDensity, by = "CHROM") %>% select(Name, win_start, win_end, isUnique_sum)

circos.genomicTrack(havePrimer, ylim = c(0, max(havePrimer$havePrimer_sum)), bg.border = NA,
                    panel.fun = function(region, value, ...) {
                      circos.genomicLines(region, value, area = TRUE, col = "red", border = NA)
                    })
circos.genomicTrack(isUnique, ylim = c(0, max(havePrimer$havePrimer_sum)), bg.border = NA,
                    panel.fun = function(region, value, ...) {
                      circos.genomicLines(region, value, area = TRUE, col = "blue", border = NA)
                    })
circos.clear()
legend(position, inset=.05, c("All primer", "Specific primer"), pch=15, col=c("red", "blue"), pt.cex = pt.cex, cex = cex, box.col = NA)
dev.off()

pdf(file = str_replace(argv[2], "primer.txt", "Primer_Distribution.pdf"), width = width, height = height)
par(mar = mar + 0.1)
#par(mar = c(5, 4, 4, 2) + 0.1)
circos.par("start.degree" = 90, track.height = track.height, track.margin = c(0, 0), cell.padding = c(0, 1.00, 0.02, 1.00))
cat(date(), ", initialization ...\n", sep = "")
circos.genomicInitialize(chromInfo %>% select(-CHROM), plotType = c("axis", "labels"))

havePrimer <- chr %>% left_join(primerDensity, by = "CHROM") %>% select(Name, win_start, win_end, havePrimer_sum)
isUnique <- chr %>% left_join(primerDensity, by = "CHROM") %>% select(Name, win_start, win_end, isUnique_sum)

circos.genomicTrack(havePrimer, ylim = c(0, max(havePrimer$havePrimer_sum)), bg.border = NA,
                    panel.fun = function(region, value, ...) {
                      circos.genomicLines(region, value, area = TRUE, col = "red", border = NA)
                    })
circos.genomicTrack(isUnique, ylim = c(0, max(havePrimer$havePrimer_sum)), bg.border = NA,
                    panel.fun = function(region, value, ...) {
                      circos.genomicLines(region, value, area = TRUE, col = "blue", border = NA)
                    })
circos.clear()
legend(position, inset=.05, c("All primer", "Specific primer"), pch=15, col=c("red", "blue"), pt.cex = pt.cex, cex = cex, box.col = NA)
dev.off()

