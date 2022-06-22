#!/bin/env Rscript

popType <- "RIL"
bulkSizeH <- 30
bulkSizeL <- 30

minDepth <- 5
maxDepth <- 200

repN <- 10000

#set.seed(123)

dltIndexMatrix_CI99up <- matrix(nrow = maxDepth, ncol = maxDepth)
dltIndexMatrix_CI99down <- matrix(nrow = maxDepth, ncol = maxDepth)
dltIndexMatrix_CI95up <- matrix(nrow = maxDepth, ncol = maxDepth)
dltIndexMatrix_CI95down <- matrix(nrow = maxDepth, ncol = maxDepth)

cat(date(), ", program start runing ...\n", sep = "")
for (depthH in minDepth:maxDepth) {
  for (depthL in minDepth:maxDepth) {
    indexH <- vector(length = repN)
    indexL <- vector(length = repN)
    dltIndex <- vector(length = repN)
    for (i in 1:repN) {
      if (popType == "RIL") {
        PH <- mean(sample(c(0,1), bulkSizeH, replace = TRUE))
        PL <- mean(sample(c(0,1), bulkSizeL, replace = TRUE))
      }else if (popType == "F2") {
        PH <- mean(sample(c(0, 0.5, 1), bulkSizeH, replace = TRUE, prob = c(0.25, 0.5, 0.25)))
        PL <- mean(sample(c(0, 0.5, 1), bulkSizeL, replace = TRUE, prob = c(0.25, 0.5, 0.25)))
      }
      indexH[i] <- mean(sample(c(0, 1), depthH, replace = TRUE, prob = c(PH, 1-PH)))
      indexL[i] <- mean(sample(c(0, 1), depthL, replace = TRUE, prob = c(PL, 1-PL)))
      dltIndex[i] <- indexL[i] - indexH[i]
    }
    dltIndexMatrix_CI95up[depthH, depthL] <- quantile(dltIndex, 0.975)
    dltIndexMatrix_CI95down[depthH, depthL] <- quantile(dltIndex, 0.025)
    dltIndexMatrix_CI99up[depthH, depthL] <- quantile(dltIndex, 0.995)
    dltIndexMatrix_CI99down[depthH, depthL] <- quantile(dltIndex, 0.005)
    #cat(depthH, depthL, 
    #    dltIndexMatrix_CI95up[depthH, depthL], dltIndexMatrix_CI95down[depthH, depthL], 
    #    dltIndexMatrix_CI99up[depthH, depthL], dltIndexMatrix_CI99down[depthH, depthL], 
    #    "\n", sep = "\t")
  }
  cat(date(), ", ", paste(depthH, maxDepth, sep = "/"), " have been done ...\n", sep = "")
}

save(dltIndexMatrix_CI95up, dltIndexMatrix_CI95down, dltIndexMatrix_CI99up, dltIndexMatrix_CI99down, 
     file = paste("QTLseqCI", 
                  paste(paste("indH", bulkSizeH, sep = ""), paste("indL", bulkSizeL, sep = ""), 
                        popType, sep = "_"), 
                  paste("Depth", minDepth, maxDepth, sep = "_"), 
                  paste("Rep", repN, sep = "_"),
                  "RData", sep = "."))
