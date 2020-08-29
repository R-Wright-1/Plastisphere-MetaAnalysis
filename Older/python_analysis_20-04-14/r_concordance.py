#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 21 15:39:03 2020

@author: robynwright
"""

compute_concordance <- function(ps_fitted_list){
  conc_df <- NULL
  for(i in 1:length(ps_fitted_list)){ # i in 1:100 comparisons
    # pval extraction
    P_df1 <- llply(.data = ps_fitted_list[[i]]$Subset1,.fun = function(method){
      if("pValMat" %in% names(method)){
        out <- method$pValMat[!is.na(method$pValMat[,1]) & method$pValMat[,1]<1,1]
      } else {
        if (class(method) == "data.frame"){
          # The CATplot function will sort the differentials of songbird with decreasing = FALSE
          out <- -abs(method[,2]) # So the first feature should be the largest in rank (that's why the minus)
          names(out) <- rownames(method)
        } else { # The same for the loadings of mixMC 
          out <- -abs(method$statInfo[,2])
          names(out) <- rownames(method$statInfo)
        }
      }
      return(out)
    }) # first half data
    P_df2 <- llply(.data = ps_fitted_list[[i]]$Subset2,.fun = function(method){
      if("pValMat" %in% names(method)){
        out <- method$pValMat[!is.na(method$pValMat[,1]) & method$pValMat[,1]<1,1]
      } else {
        if (class(method) == "data.frame"){
          # The CATplot function will sort the differentials of songbird with decreasing = FALSE
          out <- -abs(method[,2]) # So the first feature should be the largest in rank (that's why the minus)
          names(out) <- rownames(method)
        } else { # The same for the loadings of mixMC 
          out <- -abs(method$statInfo[,2])
          names(out) <- rownames(method$statInfo)
        }
      }
      return(out)
    })  # second half data
    
    nmethods <- length(names(ps_fitted_list[[i]]$Subset1))
    for(j in 1:nmethods){ # j in method names
      cat("Mehod",names(ps_fitted_list[[i]]$Subset1)[j],"with \n")
      for(k in 1:nmethods){ # k in method names again
        cat("\t",names(ps_fitted_list[[i]]$Subset1)[k],"\n")
        if(j != k){ # BMC computation
          # BMC for Subset1 data
          conc_subset1 <- data.frame(CATplot(vec1 = P_df1[[j]],vec2 = P_df1[[k]],make.plot = FALSE,maxrank = 100), 
                                     method1 = names(ps_fitted_list[[i]]$Subset1)[j], 
                                     method2 = names(ps_fitted_list[[i]]$Subset1)[k],
                                     #ndisc_0.1_method1 = length(adjP_df1[[j]]),
                                     #ndisc_0.1_method2 = length(adjP_df1[[k]]),
                                     nfeatures = ifelse(test = ("pValMat" %in% names(ps_fitted_list[[i]]$Subset1[[j]])),
                                                        yes = nrow(ps_fitted_list[[i]]$Subset1[[j]]$pValMat),
                                                        no = ifelse(test = class(ps_fitted_list[[i]]$Subset1[[j]]) == "data.frame",
                                                                    yes = nrow(ps_fitted_list[[i]]$Subset1[[j]]),
                                                                    no = nrow(ps_fitted_list[[i]]$Subset1[[j]]$statInfo))),
                                     comparison = i,
                                     subset = "1")
          # BMC for Subset2 data
          conc_subset2 <- data.frame(CATplot(vec1 = P_df2[[j]],vec2 = P_df2[[k]],make.plot = FALSE,maxrank = 100), 
                                     method1 = names(ps_fitted_list[[i]]$Subset2)[j], 
                                     method2 = names(ps_fitted_list[[i]]$Subset2)[k],
                                     #ndisc_0.1_method1 = length(adjP_df2[[j]]),
                                     #ndisc_0.1_method2 = length(adjP_df2[[k]]),
                                     nfeatures = ifelse(test = ("pValMat" %in% names(ps_fitted_list[[i]]$Subset2[[j]])),
                                                        yes = nrow(ps_fitted_list[[i]]$Subset2[[j]]$pValMat),
                                                        no = ifelse(test = class(ps_fitted_list[[i]]$Subset2[[j]]) == "data.frame",
                                                                    yes = nrow(ps_fitted_list[[i]]$Subset2[[j]]),
                                                                    no = nrow(ps_fitted_list[[i]]$Subset2[[j]]$statInfo))),
                                     comparison = i,
                                     subset = "2")
          conc <- rbind(conc_subset1,conc_subset2)
        } else {
          # WMC computed between Subset1 and Subset2
          conc <- data.frame(CATplot(vec1 = P_df1[[j]],vec2 = P_df2[[k]],make.plot = FALSE,maxrank = 100), 
                             method1 = names(ps_fitted_list[[i]]$Subset1)[j], 
                             method2 = names(ps_fitted_list[[i]]$Subset2)[k],
                             #ndisc_0.1_method1 = length(adjP_df1[[j]]),
                             #ndisc_0.1_method2 = length(adjP_df2[[k]]),
                             nfeatures = mean(ifelse(test = ("pValMat" %in% names(ps_fitted_list[[i]]$Subset1[[j]])),
                                                     yes = nrow(ps_fitted_list[[i]]$Subset1[[j]]$pValMat),
                                                     no = ifelse(test = class(ps_fitted_list[[i]]$Subset1[[j]]) == "data.frame",
                                                                 yes = nrow(ps_fitted_list[[i]]$Subset1[[j]]),
                                                                 no = nrow(ps_fitted_list[[i]]$Subset1[[j]]$statInfo))),
                                              ifelse(test = ("pValMat" %in% names(ps_fitted_list[[i]]$Subset2[[k]])),
                                                     yes = nrow(ps_fitted_list[[i]]$Subset2[[k]]$pValMat),
                                                     no = ifelse(test = class(ps_fitted_list[[i]]$Subset2[[k]]) == "data.frame",
                                                                 yes = nrow(ps_fitted_list[[i]]$Subset2[[k]]),
                                                                 no = nrow(ps_fitted_list[[i]]$Subset2[[k]]$statInfo)))),
                             comparison = i,
                             subset = "1vs2")
        }
        conc_df <- rbind(conc_df,conc)
      }
    }
  }
  return(conc_df)
}