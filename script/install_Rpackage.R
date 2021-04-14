#!/bin/env Rscript
install.packages(c("devtools", "argparser", "tidyverse", "cowplot", "ggsci"))
library(devtools)
devtools::install_github("bmansfeld/QTLseqr")
devtools::install_github('tavareshugo/windowscanr')
