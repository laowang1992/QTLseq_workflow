library(tidyverse)

getQTLseqCI <- function(df, outPrefix, popType, bulkSizeH, bulkSizeL, minDepth = 5, maxDepth = 150, repN = 10000){
  dltIndex_CI <- df %>% mutate(CI95upper = NA, CI95lower = NA, CI99upper = NA, CI99lower = NA)
  #dltIndex_CI <- tibble(HB.DP = rep(minDepth:maxDepth, times = maxDepth-minDepth+1), 
  #                      LB.DP = rep(minDepth:maxDepth, each = maxDepth-minDepth+1),
  #                      CI95upper = NA, CI95lower = NA, CI99upper = NA, CI99lower = NA)
  cat(date(), ", program start runing ...\n", sep = "")
  # 定义一个进度条
  width <- options()$width
  pb <- progress::progress_bar$new(
    format = 'Progress [:bar] :percent eta: :eta',
    total = nrow(dltIndex_CI), clear = FALSE, width = width
  )
  for (i in 1:nrow(dltIndex_CI)) {
    depthH <- dltIndex_CI[i, 1][[1]]
    depthL <- dltIndex_CI[i, 2][[1]]
    
    if (popType == "RIL") {
      PH <- rbinom(repN, bulkSizeH, 0.5) / bulkSizeH
      PL <- rbinom(repN, bulkSizeL, 0.5) / bulkSizeL
      #PH <- apply(rmultinom(repN, bulkSizeH, c(1, 1)) * c(1, 0) / bulkSizeH, 2, sum)
      #PL <- apply(rmultinom(repN, bulkSizeL, c(1, 1)) * c(1, 0) / bulkSizeL, 2, sum)
    } else if (popType == "F2") {
      PH <- apply(rmultinom(repN, bulkSizeH, c(1, 2, 1)) * c(1, 0.5, 0) / bulkSizeH, 2, sum)
      PL <- apply(rmultinom(repN, bulkSizeL, c(1, 2, 1)) * c(1, 0.5, 0) / bulkSizeL, 2, sum)
    }
    
    indexH <- rbinom(repN, depthH, PH) / depthH
    indexL <- rbinom(repN, depthL, PL) / depthL
    dltIndex <- indexH - indexL
    
    dltIndex_CI[i, c("CI95upper", "CI95lower", "CI99upper", "CI99lower")] <- t(
      c(quantile(dltIndex, 0.975), quantile(dltIndex, 0.025), 
        quantile(dltIndex, 0.995), quantile(dltIndex, 0.005))
    )
    
    # 打印进度条
    pb$tick()
  }
  cat(date(), ", finish, export data ...\n", sep = "")
  
  
  df <- gather(dltIndex_CI, type, CI, -HB.DP, -LB.DP) %>% 
    mutate(level = if_else(str_detect(type, "95"), "95 CI", "99 CI"),
           direction = if_else(str_detect(type, "upper"), "Upper CI", "Lower CI"))
  df$direction <- factor(df$direction, levels = c("Upper CI", "Lower CI"))
  P_CI <- ggplot(df, aes(x = HB.DP, y = LB.DP, fill = CI)) + 
    geom_tile() + 
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_distiller(palette = "RdBu") +
    labs(title = paste(outPrefix, "CI", 
                       paste(paste("indH", bulkSizeH, sep = ""), 
                             paste("indL", bulkSizeL, sep = ""), 
                             popType, sep = "_"), 
                       paste("Depth", minDepth, maxDepth, sep = "_"), 
                       paste("Rep", repN, sep = "_"), sep = "."),
         fill = NULL) +
    facet_grid(direction ~ level) +
    theme_half_open() +
    theme(strip.background = element_rect(fill = "#90EE90"),
          plot.title = element_text(hjust = 0.5))
  ggsave(P_CI, filename = paste(outPrefix, "QTLseqCI", 
                                paste(paste("indH", bulkSizeH, sep = ""), 
                                      paste("indL", bulkSizeL, sep = ""), 
                                      popType, sep = "_"), 
                                paste("Depth", minDepth, maxDepth, sep = "_"), 
                                paste("Rep", repN, sep = "_"),
                                "pdf", sep = "."),
         width = 6.5, height = 6)
  ggsave(P_CI, filename = paste(outPrefix, "QTLseqCI", 
                                paste(paste("indH", bulkSizeH, sep = ""), 
                                      paste("indL", bulkSizeL, sep = ""), 
                                      popType, sep = "_"), 
                                paste("Depth", minDepth, maxDepth, sep = "_"), 
                                paste("Rep", repN, sep = "_"),
                                "png", sep = "."),
         width = 6.5, height = 6, units = "in", dpi = 500)
  #save(dltIndex_CI, 
  #     file = paste("QTLseqCI", 
  #                  paste(paste("indH", bulkSizeH, sep = ""), 
  #                        paste("indL", bulkSizeL, sep = ""), 
  #                        popType, sep = "_"), 
  #                  paste("Depth", minDepth, maxDepth, sep = "_"), 
  #                  paste("Rep", repN, sep = "_"),
  #                  "RData", sep = "."))
  return(dltIndex_CI)
}


getQTL <- function(data, chr = "CHROM", pos = "POS", 
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
  
  # Peak information
  interval <- interval %>% mutate(PeakPos = numeric(nrow(interval)), PeakDeltaIndex = numeric(nrow(interval)))
  for (i in seq_along(interval$CHROM)) {
    peak <- df %>% filter(chr == interval$CHROM[i], pos >= interval$Start[i], pos <= interval$End[i]) %>% arrange(desc(index))
    if (peak$index[1] > peak$CIupper[1]) {
      interval$PeakPos[i] <- peak %>% head(1) %>% pull(pos)
      interval$PeakDeltaIndex[i] <- peak %>% head(1) %>% pull(index)
    }else{
      interval$PeakPos[i] <- peak %>% tail(1) %>% pull(pos)
      interval$PeakDeltaIndex[i] <- peak %>% tail(1) %>% pull(index)
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
  interval <- read_csv(file = paste(outPrefix, paste(CI, "CI", sep = ""), "csv", sep = "."), show_col_types = FALSE)
  if (nrow(interval) > 0) {
    for (i in unique(interval$CHROM)) {
      d <- df %>% filter(CHROM == i, nSNPs >= minN) %>% 
        dplyr::select(CHROM, POS = POS, CIupper = all_of(paste("CI", CI, "upper", sep = "")), 
                      CIupper = all_of(paste("CI", CI, "upper", sep = "")), 
                      CIlower = all_of(paste("CI", CI, "lower", sep = "")), 
                      delta.index, ED, ED4)
      inter <- interval %>% filter(CHROM == i) %>% left_join(chr, by = "CHROM")
      chrLen <- len$Len[len$CHROM == i]
      p_index <- ggplot() +
        geom_rect(data = inter, aes(xmin = Start, xmax = End, ymin = -Inf, ymax = Inf),fill = '#FF3300', color = "#FF3300") +
        geom_line(data = d, aes(x = POS, y = CIupper), color = "gray60") + 
        geom_line(data = d, aes(x = POS, y = CIlower), color = "gray60") +
        geom_line(data = d, aes(x = POS, y = delta.index), color = "blue") +
        scale_x_continuous(limits = c(0, chrLen)) +
        scale_y_continuous(limits = c(-1, 1)) +
        geom_hline(aes(yintercept=0)) + 
        labs(x = chr$LABEL[chr$CHROM == i], y = "Delta SNP index") +
        theme_half_open()
      ggsave(p_index, filename = paste(outPrefix, chr$LABEL[chr$CHROM == i], paste(CI, "CI", sep = ""), "png", sep = "."), height = 3.5, width = 4.5, dpi = 500)
      ggsave(p_index, filename = paste(outPrefix, chr$LABEL[chr$CHROM == i], paste(CI, "CI", sep = ""), "pdf", sep = "."), height = 3.5, width = 4.5)
      
      p_ed <- ggplot(d, aes(x = POS, y = ED)) + 
        geom_line(color = "blue") + 
        scale_x_continuous(limits = c(0, chrLen)) +
        scale_y_continuous(limits = c(0, max(d$ED, na.rm = TRUE)*1.05)) +
        labs(x = chr$LABEL[chr$CHROM == i], y = "ED") +
        theme_half_open()
      ggsave(p_ed, filename = paste(outPrefix, chr$LABEL[chr$CHROM == i], "ED", "png", sep = "."), height = 3.5, width = 4.5, dpi = 500)
      ggsave(p_ed, filename = paste(outPrefix, chr$LABEL[chr$CHROM == i], "ED", "pdf", sep = "."), height = 3.5, width = 4.5)
      
      p_ed4 <- ggplot(d, aes(x = POS, y = ED4)) + 
        geom_line(color = "blue") + 
        scale_x_continuous(limits = c(0, chrLen)) +
        scale_y_continuous(limits = c(0, max(d$ED4, na.rm = TRUE))) +
        labs(x = chr$LABEL[chr$CHROM == i], y = bquote(ED^4)) +
        theme_half_open()
      ggsave(p_ed4, filename = paste(outPrefix, chr$LABEL[chr$CHROM == i], "ED4", "png", sep = "."), height = 3.5, width = 4.5, dpi = 500)
      ggsave(p_ed4, filename = paste(outPrefix, chr$LABEL[chr$CHROM == i], "ED4", "pdf", sep = "."), height = 3.5, width = 4.5)
      
      print(paste(chr$LABEL[chr$CHROM == i], "has been done...", sep = " "))
    }
  }
}

# A sliding window function in R.
# Author: Wang Pengfei <wangpf0608@126.com>

slidingWindow <- function(df, winSize, winStep, groups, position, values, fun){
  # winStep应小于等于winSize
  if (winStep>winSize) {
    stop("winStep should not be bigger than winSize!")
  }
  df <- df %>% dplyr::select(groups = all_of(groups), position = all_of(position), all_of(values))
  chrLen <- df %>% group_by(groups) %>% summarise(Len = max(position))
  # 生成windows
  wid <- tibble()
  for (i in 1:nrow(chrLen)) {
    if (chrLen$Len[i]<=winSize) {
      start <- 1
      end <- winSize
    }else{
      start = seq(from = 1, to = chrLen$Len[i]-winSize+winStep+1, by = winStep)
      end = seq(from = winSize, to = chrLen$Len[i]+winStep, by = winStep)
    }
    wid_tmp <- tibble(groups = chrLen$groups[i], Start = start, End = end)
    wid <- rbind(wid, wid_tmp)
  }
  # 添加统计列（N，sum）
  # 分窗口统计
  cat(date(), ", Sliding window statistics start ...\n", sep = "")
  # 定义一个进度条
  width <- options()$width
  pb <- progress::progress_bar$new(
    format = 'Progress [:bar] :percent eta: :eta',
    total = nrow(wid), clear = FALSE, width = width
  )
  for (i in 1:nrow(wid)) {
    df_tmp <- df %>% filter(groups == wid$groups[i], position >= wid$Start[i], position <= wid$End[i]) %>% na.omit()
    wid[i, "N"] <- nrow(df_tmp)
    x <- apply(df_tmp[, values], 2, fun)
    wid[i, paste(values, fun, sep = "_")] <- t(x[values])
    
    # 打印一个进度条
    ## 使用progress扩展包打印进度条
    pb$tick()
    ## 以下是一个手动打印进度条的代码
    #########################################################
    #cat('[', paste0(rep('#', i/nrow(wid)*width), collapse=''),
    #    paste0(rep('-', width - i/nrow(wid)*width), collapse=''),
    #    ']',
    #    round(i/nrow(wid)*100),'%')
    #if(i==nrow(wid))cat('\n', date(), ' DONE!\n', sep = "")
    #else cat('\r')
    ##########################################################
  }
  cat(date(), ", Done!\n", sep = "")
  colnames(wid)[1:3] <- c(groups, "win_start", "win_end")
  return(wid)
}

addUp <- function(df, len = NULL, group, pos, band = 0.01){
  df_tmp <- df
  df <- df %>% ungroup() %>% dplyr::select(group = all_of(group), all_of(pos))
  # 取最大值
  if (is.null(len)) {
    len <- df %>% mutate(max = apply(df[,-1], 1, max)) %>% group_by(group) %>% summarise(Len = max(max))
  }else{
    len <- len %>% dplyr::select(group = all_of(group), Len)
  }
  accu <- rep(0, rep(nrow(len)))
  for (i in seq_along(accu)[-1]) {
    accu[i] <- sum(len$Len[1:i-1]) + sum(len$Len)*band*(i-1)
  }
  names(accu) <- len$group
  breaks <- accu + len$Len/2
  gaps <- accu[-1] - sum(len$Len)*band/2
  labels <- names(accu)
  
  for (p in pos) {
    dd <- df %>% select(group, p = all_of(p)) %>% mutate(p_addUp = 0)
    for (c in names(accu)) {
      dd <- dd %>% mutate(p_addUp = if_else(group == c, p+accu[c], p_addUp))
    }
    colnames(dd) <- c("group", p, paste(p, "addUp", sep = "_"))
    df <- df %>% left_join(dd, by = c("group", p))
  }
  
  colnames(df)[1] <- group
  df <- df %>% left_join(df_tmp, by = c(group, pos))
  outList <- list(df = df, breaks = breaks, labels = labels, gaps = gaps)
  return(outList)
}
