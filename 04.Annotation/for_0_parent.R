selectSample_and_filterGT <- function(data, highP, lowP, highB, lowB, minGQ){
  #
  dd1 <- data %>%
    select(CHROM, POS, REF, ALT,
           all_of(c(paste(highB, "GT", sep = "."), paste(highB, "AD", sep = "."), paste(highB, "GQ", sep = "."))),
           all_of(c(paste(lowB, "GT", sep = "."), paste(lowB, "AD", sep = "."), paste(lowB, "GQ", sep = ".")))) %>%
    na.omit() %>% rename(highBulk.GT = paste(highB, "GT", sep = "."),
                         lowBulk.GT = paste(lowB, "GT", sep = "."),
                         highBulk.AD = paste(highB, "AD", sep = "."),
                         lowBulk.AD = paste(lowB, "AD", sep = "."),
                         highBulk.GQ = paste(highB, "GQ", sep = "."),
                         lowBulk.GQ = paste(lowB, "GQ", sep = "."))

  # 去除无多态性位点
  dd2 <- dd1 %>% filter(!((highBulk.GT == paste(REF, REF, sep = "/") & lowBulk.GT == paste(REF, REF, sep = "/")) |
                            (highBulk.GT == paste(ALT, ALT, sep = "/") & lowBulk.GT == paste(ALT, ALT, sep = "/"))))

  # 将highBulk和lowBulk的AD分开，并计算DP值
  dd3 <- dd2 %>%
    separate(highBulk.AD, c("highBulk.AD_0", "highBulk.AD_1"), sep = ",", convert = TRUE) %>%
    separate(lowBulk.AD, c("lowBulk.AD_0", "lowBulk.AD_1"), sep = ",", convert = TRUE) %>%
    mutate(highBulk.DP = highBulk.AD_0 + highBulk.AD_1,
           lowBulk.DP = lowBulk.AD_0 + lowBulk.AD_1) # 因为highBulk.AD_0 + highBulk.AD_1有时并不等于原来的highBulk.DP并且会导致后面代码报错而且很难发现原因，因此这里重新计算覆盖原来的highBulk.DP

  # 过滤HB index和LB index同时高于某值或同时低于某值的位点
  # 这里的数值可以改成一个参数，后续修改可以考虑
  dd4 <- dd3 %>% filter(!(((highBulk.AD_0 / highBulk.DP < 0.3) &
                             (lowBulk.AD_0 / lowBulk.DP < 0.3)) |
                            ((highBulk.AD_0 / highBulk.DP > 0.7) &
                               (lowBulk.AD_0 / lowBulk.DP > 0.7))))

  # 根据亲本、混池基因型，将混池的AD_0和AD_1重新分配
  dd5 <- dd4 %>%
    dplyr::select(CHROM, POS, REF, ALT,
                  HB.REF.AD = highBulk.AD_0, HB.ALT.AD = highBulk.AD_1, HB.DP = highBulk.DP,
                  LB.REF.AD = lowBulk.AD_0, LB.ALT.AD = lowBulk.AD_1, LB.DP = lowBulk.DP)
  return(dd5)
}

depth_statistics <- function(data, outPrefix, highP, lowP, highB, lowB){
  # 统计个样本的DP分布密度
  dp <- data %>% dplyr::select(HB = HB.DP, LB = LB.DP) %>%
    gather(key = "sample", value = "depth")
  dp$sample <- factor(dp$sample, levels = c("HB", "LB"), labels = c(highB, lowB))
  upper <- dp %>% group_by(sample) %>% summarise(ave = mean(depth), sd = sd(depth)) %>% mutate(upper = ave+6*sd) %>% pull(upper) %>% max()
  P_dp <- ggplot(dp, aes(x = depth)) +
    geom_histogram(aes(y = after_stat(density), fill = sample), binwidth = 2) +
    geom_density() +
    scale_x_continuous(limits = c(0, upper)) +
    theme_half_open() +
    facet_wrap(~sample, nrow = 1)
  ggsave(P_dp, filename = paste(outPrefix, "depth_density.pdf", sep = "."), height = 3.5, width = 6)
  ggsave(P_dp, filename = paste(outPrefix, "depth_density.png", sep = "."), height = 3.5, width = 6, dpi = 500)
}

filterDP <- function(data, minHPdp, minLPdp, minHBdp, minLBdp){
  # 根据深度过滤，上限由以前的自定义改为设置成ave+3*sd
  dp <- data %>% dplyr::select(HB = HB.DP, LB = LB.DP) %>%
    gather(key = "sample", value = "depth")
  stat <- dp %>% group_by(sample) %>%
    summarise(ave = mean(depth), sd = sd(depth), upper = ave+3*sd, outer = sum(depth>upper), `outerRate` = paste(round(x = outer/n()*100, digits = 3), "%"))
  print(stat)
  maxHBdp <<- ceiling(stat$upper[stat$sample %in% "HB"])
  maxLBdp <<- ceiling(stat$upper[stat$sample %in% "LB"])
  # 根据深度过滤
  df <- data %>% filter(HB.DP > minHBdp, HB.DP < maxHBdp,
                        LB.DP > minLBdp, LB.DP < maxLBdp)
  df
}

SNP_distribution <- function(data, outPrefix, chr){
  # 计算每条染色体位点个数
  SNPnumber <- data %>% group_by(CHROM) %>% count()
  write_tsv(x = SNPnumber, file = paste(outPrefix, "SNP_number_per_chr.txt", sep = "."))
  write_csv(x = SNPnumber, file = paste(outPrefix, "SNP_number_per_chr.csv", sep = "."))

  # SNP分布图
  options(scipen = 200)
  colourCount = nrow(chr)
  getPalette = colorRampPalette(brewer.pal(8, "Set1"))
  chr$LABEL <- factor(chr$LABEL, levels = chr$LABEL)
  Phist <- chr %>%
    left_join(data, by = "CHROM") %>%
    ggplot(aes(x = POS)) +
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
}

export_dp <- function(data, outPrefix, highP, lowP, highB, lowB){
  # export allele Depth data
  outtb <- data %>%
    dplyr::select(CHROM, POS, REF, ALT, HB.REF.AD, HB.ALT.AD, LB.REF.AD, LB.ALT.AD)
  colnames(outtb) <- c("CHROM", "POS", "REF", "ALT",
                       paste(highB, "ref", "depth", sep = "."),
                       paste(highB, "alt", "depth", sep = "."),
                       paste(lowB, "ref", "depth", sep = "."),
                       paste(lowB, "alt", "depth", sep = "."))
  write_tsv(x = outtb, file = paste(outPrefix, "Depth_information.txt", sep = "."))
  write_csv(x = outtb, file = paste(outPrefix, "Depth_information.csv", sep = "."))
}

calc_index_etc <- function(data, CI = NULL, outPrefix, minHBdp, minLBdp, popType, bulkSizeH, bulkSizeL, winSize, winStep){
  df <- data %>%
    dplyr::mutate(HB.index = HB.REF.AD / HB.DP,
                  LB.index = LB.REF.AD / LB.DP,
                  ED = sqrt((HB.index - LB.index)^2 + ((1-HB.index) - (1-LB.index))^2))
  # 滑窗统计，windowscanr这个包好久不更新，说不定那天就不弄用了，还是自己写一个滑窗统计函数吧
  slidwin <- slidingWindow(df = df,
                           winSize = winSize,
                           winStep = winStep,
                           groups = "CHROM",
                           position = "POS",
                           values = c("ED"),
                           fun = "mean") %>% as_tibble() %>%
    dplyr::select(CHROM, win_start, win_end,
                  ED= ED_mean, nSNPs = N) %>%
    dplyr::mutate(ED4 = ED^4, POS = win_start/2 + win_end/2)
  write_tsv(x = slidwin %>% dplyr::select(-POS), file = paste(outPrefix, "SlidingWindow.txt", sep = "."))
  write_csv(x = slidwin %>% dplyr::select(-POS), file = paste(outPrefix, "SlidingWindow.csv", sep = "."))
  slidwin
}

export_figure <- function(data, outPrefix, chr, len, minN, highB, lowB, width, height,
                          color = c("#4197d8","#f8c120", "#413496", "#495226", "#d60b6f", "#e66519", "#d581b7", "#83d3ad", "#7c162c", "#26755d")){
  COLOR <- chr %>% mutate(color = rep(color, len = nrow(chr)))
  newLen <- chr %>% left_join(len, by = "CHROM") %>% select(LABEL, Len)

  dataforPlot <- data %>% filter(nSNPs >= minN) %>% right_join(chr, by = "CHROM") %>%
    addUp(len = newLen, group = "LABEL", pos = "POS", band = 0.005)


  p <- ggplot(dataforPlot$df, aes(x = POS_addUp, group = CHROM, color = CHROM)) +
    geom_vline(xintercept = dataforPlot$gaps, linetype = "dashed", color = "gray") +
    scale_x_continuous(breaks = dataforPlot$breaks, labels = dataforPlot$labels, expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_color_manual(breaks = COLOR$CHROM, values = COLOR$color) +
    theme_half_open() +
    theme(legend.position = "NULL",
          #axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title.x = element_blank())

  #
  p_ED <- p +
    coord_cartesian(ylim = c(0, max(dataforPlot$df$ED, na.rm = TRUE)*1.05)) +
    labs(y = "ED")
  p_ED_L <- p_ED + geom_line(aes(y = ED), linewidth = 1)
  ggsave(p_ED_L, filename = paste(outPrefix, "ED.line.pdf", sep = "."), width = width, height = height)
  ggsave(p_ED_L, filename = paste(outPrefix, "ED.line.png", sep = "."), width = width, height = height, dpi = 500)
  p_ED_P <- p_ED + geom_point(aes(y = ED), size = 1.5)
  ggsave(p_ED_P, filename = paste(outPrefix, "ED.point.pdf", sep = "."), width = width, height = height)
  ggsave(p_ED_P, filename = paste(outPrefix, "ED.point.png", sep = "."), width = width, height = height, dpi = 500)

  #
  p_ED4 <- p +
    coord_cartesian(ylim = c(0-max(dataforPlot$df$ED4, na.rm = TRUE)*0.015, max(dataforPlot$df$ED4, na.rm = TRUE)*1.05)) +
    labs(y = bquote(ED^4))
  p_ED4_L <- p_ED4 + geom_line(aes(y = ED4), linewidth = 1)
  ggsave(p_ED4_L, filename = paste(outPrefix, "ED4.line.pdf", sep = "."), width = width, height = height)
  ggsave(p_ED4_L, filename = paste(outPrefix, "ED4.line.png", sep = "."), width = width, height = height, dpi = 500)
  p_ED4_P <- p_ED4 + geom_point(aes(y = ED4), size = 1.5)
  ggsave(p_ED4_P, filename = paste(outPrefix, "ED4.point.pdf", sep = "."), width = width, height = height)
  ggsave(p_ED4_P, filename = paste(outPrefix, "ED4.point.png", sep = "."), width = width, height = height, dpi = 500)
}


