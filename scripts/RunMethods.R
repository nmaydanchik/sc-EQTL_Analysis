RunMethods <- function(mr.obj, allMethods = FALSE) {
  
  if (allMethods) est = matrix(nrow = 10, ncol = 3) else est = matrix(nrow = 7, ncol = 3)
  colnames(est) = c("Method", "Estimate", "p-value")
  
  complete = 0
  
  # 1. IVW
  res = mr_ivw(mr.obj)
  est[complete + 1, 1] = "IVW"
  est[complete + 1, 2] = res$Estimate
  est[complete + 1, 3] = res$Pvalue
  complete = complete + 1
  print(paste(complete, "Methods completed"))
  rm(res)
  
  # 2. median
  res = mr_median(mr.obj)
  est[complete + 1, 1] = "Median"
  est[complete + 1, 2] = res$Estimate
  est[complete + 1, 3] = res$Pvalue
  complete = complete + 1
  print(paste(complete, "Methods completed"))
  rm(res)
  
  # 3. mode
  res = mr_mbe(mr.obj)
  est[complete + 1, 1] = "Mode"
  est[complete + 1, 2] = res$Estimate
  est[complete + 1, 3] = res$Pvalue
  complete = complete + 1
  print(paste(complete, "Methods completed"))
  rm(res)
  
  # 4. robust
  res = mr_ivw(mr.obj,"random", robust = TRUE)
  est[complete + 1, 1] = "Robust"
  est[complete + 1, 2] = res$Estimate
  est[complete + 1, 3] = res$Pvalue
  complete = complete + 1
  print(paste(complete, "Methods completed"))
  rm(res)
  
  # 5 Egger
  res = mr_egger(mr.obj)
  est[complete + 1, 1] = "Egger"
  est[complete + 1, 2] = res$Estimate
  est[complete + 1, 3] = res$Pvalue.Est
  complete = complete + 1
  print(paste(complete, "Methods completed"))
  rm(res)
  
  # 6. contamination mixture
  res = mr_conmix(mr.obj)
  est[complete + 1, 1] = "ConMix"
  est[complete + 1, 2] = res$Estimate
  est[complete + 1, 3] = res$Pvalue
  complete = complete + 1
  print(paste(complete, "Methods completed"))
  rm(res)
  
  # 7. MR-RAPS
  res = mr.raps(mr.obj$betaX, mr.obj$betaY, mr.obj$betaXse, mr.obj$betaYse)
  est[complete + 1, 1] = "RAPS"
  est[complete + 1, 2] = res$beta.hat
  est[complete + 1, 3] = res$beta.p.value
  complete = complete + 1
  print(paste(complete, "Methods completed"))
  rm(res)
  
  if (allMethods) {
    # 8. MR-Grapple
    grapple.df = data.frame(SNP = mr.obj$snps, gamma_exp1 = mr.obj$betaX, gamma_out1 = mr.obj$betaY, se_exp1 = mr.obj$betaXse, se_out1 = mr.obj$betaYse)
    res = grappleRobustEst(grapple.df, plot.it = F)
    est[complete + 1, 1] = "Grapple"
    est[complete + 1, 2] = res$beta.hat
    est[complete + 1, 3] = res$beta.p.value
    complete = complete + 1
    print(paste(complete, "Methods completed"))
    rm(res, grapple.df)    
    
    # 9. MR-Mix
    res = MRMix(mr.obj$betaX, mr.obj$betaY, mr.obj$betaXse, mr.obj$betaYse)
    est[complete + 1, 1] = "Mix"
    est[complete + 1, 2] = res$theta
    est[complete + 1, 3] = res$pvalue_theta
    complete = complete + 1
    print(paste(complete, "Methods completed"))
    rm(res)
    
    
    # 10. MR-Lasso
    res = mr_lasso(mr.obj)
    est[complete + 1, 1] = "Lasso"
    est[complete + 1, 2] = res$Estimate
    est[complete + 1, 3] = res$Pvalue
    complete = complete + 1
    print(paste(complete, "Methods completed"))
    rm(res)
  }
  
  return(est)
}
