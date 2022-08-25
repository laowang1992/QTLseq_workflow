library(tidyverse)
getQTL <- function(data, chr = "CHROM", pos = "win_mid", 
                   index = "delta.index", nSNPs = "nSNPs", CI = 95, n = 10, 
                   export = TRUE, filename = "SignificantQTL.csv"){
  ## 整理格式
  df <- data %>% 
    dplyr::select(chr = all_of(chr), pos = all_of(pos), index = all_of(index), nSNPs = all_of(nSNPs),
                  CIupper = all_of(paste("CI", CI, "upper", sep = "")),
                  CIlower = all_of(paste("CI", CI, "lower", sep = "")))
  ## 初始化结果变量
  interval <- tibble(CHROM = character(0), Start = numeric(0), End = numeric(0), Length = numeric(0))
  ##
  for (chromosome in unique(df$chr)) {
    # 提取子集
    subdf <- df %>% filter(chr == chromosome & nSNPs >= n)
    #subdf <- slidwin %>% filter(CHROM == "scaffoldA03")
    # 初始化变量
    start <- NULL
    end <- NULL
    sig <- FALSE
    for (i in seq_along(subdf$chr)) {
      if ((subdf$index[i]>subdf$CIupper[i] | subdf$index[i]<subdf$CIlower[i]) & (!sig)) {
        start <- subdf$pos[i]
        sig <- TRUE
      }
      if (((subdf$index[i]<subdf$CIupper[i] & subdf$index[i]>subdf$CIlower[i]) | i == nrow(subdf)) & sig) {
        if (i == nrow(subdf)) {
          end = subdf$pos[i]
        } else {
          end = subdf$pos[i-1]
        }
        sig = FALSE
        new_interval <- tibble(CHROM = chromosome, Start = start, End = end, Length = end - start)
        interval <- rbind(interval, new_interval)
      }
    }
  }
  if (export == TRUE) {
    write_csv(x = interval, file = filename)
  }
  return(interval)
}



plotIndex <- function(df, chr, len, CI=NA, nSNPs=NA, n=5, index, band = 0.02, 
                      color=c("#4197d8","#f8c120", "#413496", "#495226", 
                              "#d60b6f", "#e66519", "#d581b7", "#83d3ad", "#7c162c", "#26755d"), 
                      ylim = NULL, ylab="Delta SNP index", size=1, axis.size=1, 
                      axis.lab.size=1, axis.title.size=1, type="p"){
  chr <- chr %>% mutate(COLOR = rep(color, len = nrow(chr)))
  len <- chr %>% left_join(len, by = "CHROM") %>% dplyr::select(CHROM, Len)
  addUp <- vector(mode = "integer", length = nrow(len))
  s <- 0
  for (i in seq_along(addUp)) {
    addUp[i] <- s
    s = s + len$Len[i] + band*sum(len$Len)
  }
  names(addUp) <- len$CHROM
  if (is.na(CI)) {
    if (type == "p") {
      cat("here, no CI and p\n")
      df <- chr %>% left_join(df, by = "CHROM") %>% left_join(len, by = "CHROM") %>%
        dplyr::select(CHROM, LABEL, COLOR, POS, "index"=all_of(index), "nSNPs"=all_of(nSNPs))
      x <- rep(addUp, table(df$CHROM))
      df <- df %>% mutate(POS = POS+x) %>% filter(nSNPs >= n)
      # 坐标轴
      breaks <- addUp + len$Len/2
      #
      plot(x = df$POS, y = df$index, col = df$COLOR, type = "p", pch = 20, cex = size,
           ylim = ylim, axes=FALSE, xlab = NA, ylab = ylab, font.lab = 2, 
           cex.lab = axis.title.size)
      axis(side = 1, at = breaks, labels = chr$LABEL, lwd = axis.size, lwd.ticks = axis.size, 
           font = 2, cex.axis = axis.lab.size)
      axis(side = 2, las="2", lwd = axis.size, lwd.ticks = axis.size, font = 2,
           cex.axis = axis.lab.size)
    } else if (type == "l") {
      cat("here, no CI and l\n")
      df <- chr %>% left_join(df, by = "CHROM") %>% left_join(len, by = "CHROM") %>%
        dplyr::select(CHROM, LABEL, COLOR, POS, "index"=all_of(index), "nSNPs"=all_of(nSNPs))
      x <- rep(addUp, table(df$CHROM))
      df <- df %>% mutate(POS = POS+x) %>% filter(nSNPs >= n)
      # 坐标轴
      breaks <- addUp + len$Len/2
      #
      plot(0, 0, xlim = c(0, sum(len$Len)+band*sum(len$Len)*(nrow(len)-1)), ylim = ylim, type = "n", 
           axes = F, xlab = NA, ylab = ylab, font.lab = 2, cex.lab = axis.title.size)
      axis(side = 1, at = breaks, labels = chr$LABEL, lwd = axis.size, lwd.ticks = axis.size, 
           font = 2, cex.axis = axis.lab.size)
      axis(side = 2, las="2", lwd = axis.size, lwd.ticks = axis.size, font = 2,
           cex.axis = axis.lab.size)
      for (i in chr$CHROM) {
        d <- df %>% filter(CHROM == i)
        lines(x = d$POS, y = d$index, col = d$COLOR, lwd = size)
      }
    }
  } else {
    if (type == "p") {
      cat("here, CI and p\n")
      df <- chr %>% left_join(df, by = "CHROM") %>% left_join(len, by = "CHROM") %>%
        dplyr::select(CHROM, LABEL, COLOR, POS, "index"=all_of(index), "nSNPs"=all_of(nSNPs), 
                      "CIupper" = all_of(paste(CI, "upper", sep = "")),
                      "CIlower" = all_of(paste(CI, "lower", sep = "")))
      x <- rep(addUp, table(df$CHROM))
      df <- df %>% mutate(POS = POS+x) %>% filter(nSNPs >= n)
      # 坐标轴
      breaks <- addUp + len$Len/2
      #
      plot(x = df$POS, y = df$index, col = df$COLOR, type = "p", pch = 20, cex = size,
           ylim = ylim, axes=FALSE, xlab = NA, ylab = ylab, font.lab = 2, cex.lab = axis.title.size)
      axis(side = 1, at = breaks, labels = chr$LABEL, lwd = axis.size, lwd.ticks = axis.size, 
           font = 2, cex.axis = axis.lab.size)
      axis(side = 2, las="2", lwd = axis.size, lwd.ticks = axis.size, font = 2,
           cex.axis = axis.lab.size)
      for (i in chr$CHROM) {
        d <- df %>% filter(CHROM == i)
        lines(x = d$POS, y = d$CIupper, col = "gray")
        lines(x = d$POS, y = d$CIlower, col = "gray")
      }
    } else if (type == "l") {
      cat("here, CI and l\n")
      df <- chr %>% left_join(df, by = "CHROM") %>% left_join(len, by = "CHROM") %>%
        dplyr::select(CHROM, LABEL, COLOR, POS, "index"=all_of(index), "nSNPs"=all_of(nSNPs), 
                      "CIupper" = all_of(paste(CI, "upper", sep = "")),
                      "CIlower" = all_of(paste(CI, "lower", sep = "")))
      x <- rep(addUp, table(df$CHROM))
      df <- df %>% mutate(POS = POS+x) %>% filter(nSNPs >= n)
      # 坐标轴
      breaks <- addUp + len$Len/2
      #
      plot(0, 0, xlim = c(0, sum(len$Len)+band*sum(len$Len)*(nrow(len)-1)), ylim = ylim, type = "n", 
           axes = F, xlab = NA, ylab = ylab, font.lab = 2, cex.lab = axis.title.size)
      axis(side = 1, at = breaks, labels = chr$LABEL, lwd = axis.size, lwd.ticks = axis.size, 
           font = 2, cex.axis = axis.lab.size)
      axis(side = 2, las="2", lwd = axis.size, lwd.ticks = axis.size, font = 2,
           cex.axis = axis.lab.size)
      for (i in chr$CHROM) {
        d <- df %>% filter(CHROM == i)
        lines(x = d$POS, y = d$CIupper, col = "gray")
        lines(x = d$POS, y = d$CIlower, col = "gray")
        lines(x = d$POS, y = d$index, col = d$COLOR, lwd = size)
      }
    }
  }
}

plotTargetChrom <- function(df, CI = 95, minN = 0, outPrefix, chr, len){
  interval <- read_csv(file = paste(outPrefix, paste(CI, "CI", sep = ""), "csv", sep = "."))
  if (nrow(interval) > 0) {
    for (i in unique(interval$CHROM)) {
      d <- df %>% filter(CHROM == i, nSNPs >= minN) %>% 
        dplyr::select(CHROM, POS = win_mid, CIupper = all_of(paste("CI", CI, "upper", sep = "")), 
                      CIupper = all_of(paste("CI", CI, "upper", sep = "")), 
                      CIlower = all_of(paste("CI", CI, "lower", sep = "")), 
                      delta.index)
      inter <- interval %>% filter(CHROM == i) %>% left_join(chr, by = "CHROM")
      chrLen <- len$Len[len$CHROM == i]
      p <- ggplot() +
        geom_rect(data = inter, aes(xmin = Start, xmax = End, ymin = -Inf, ymax = Inf),fill = '#FF3300', color = "#FF3300") +
        geom_line(data = d, aes(x = POS, y = CIupper), color = "gray60") + 
        geom_line(data = d, aes(x = POS, y = CIlower), color = "gray60") +
        geom_line(data = d, aes(x = POS, y = delta.index), color = "blue") +
        scale_x_continuous(limits = c(0, chrLen)) +
        scale_y_continuous(limits = c(-1, 1) +
        geom_hline(aes(yintercept=0)) + 
        labs(x = chr$LABEL[chr$CHROM == i], y = "Delta SNP index") +
        theme_half_open()
      ggsave(p, filename = paste(outPrefix, chr$LABEL[chr$CHROM == i], paste(CI, "CI", sep = ""), "png", sep = "."), height = 3.5, width = 4.5, dpi = 500)
      ggsave(p, filename = paste(outPrefix, chr$LABEL[chr$CHROM == i], paste(CI, "CI", sep = ""), "pdf", sep = "."), height = 3.5, width = 4.5)
      print(paste(chr$LABEL[chr$CHROM == i], "has been done...", sep = " "))
    }
  }
}
