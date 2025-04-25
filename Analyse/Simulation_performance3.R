## ########################################################################## ##
###################     Simulation Study: results           ####################
##                          Author: KS                                        ##
##                      Date: October 2023                                    ##
##                        Out: AUC (ROC)                                      ##
## ########################################################################## ##


#################   Calculation of performance parameter    ####################

# Performance measures

# bias, MSE, cov, power Step 1
auc_vars <- grep("AUC\\.", names(res), value = TRUE)
cil_vars <- grep("CIL", names(res), value = FALSE)
ciu_vars <- grep("CIU", names(res), value = FALSE)

df <- data.frame(matrix(ncol = length(auc_vars)*4, nrow = nrow(res)))
names1 <- paste("bias", auc_vars, sep = ".")
names2 <- paste("RMSE", auc_vars, sep = ".")
names4 <- paste("cov", auc_vars, sep = ".")
names5 <- paste("power", auc_vars, sep = ".")

colnames(df) <- c(names1,names2,names4,names5)
nvars <- length(auc_vars)

for (j in 1:length(auc_vars)){
  df[,j] <- res[,auc_vars[j]] - res$AUC_0 
  df[,j+nvars] <- df[,j]^2
  cil <- cil_vars[j]
  ciu <- ciu_vars[j]
  df[,j+nvars*2] <- 1*(res[,cil]<=res$AUC_0 & res$AUC_0<=res[,ciu])
  df[,j+nvars*3] <- 1*(res$AUC_min<res[,cil])
}

raw_res <- cbind(res,df)

# average over nsim
scenario <- res[res$sim==1,fix_col] # dataframe for summary results
average <- data.frame(matrix(ncol = ncol(df)+length(auc_vars)*10, nrow = nrow(scenario))) # dataframe to store averages
names3 <- paste("empSE", auc_vars, sep = ".")
names6 <- paste("MCE.bias", auc_vars, sep = ".")
names7 <- paste("MCE.MSE", auc_vars, sep = ".")
names8 <- paste("MCE.empSE", auc_vars, sep = ".")
names9 <- paste("MCE.cov", auc_vars, sep = ".")
names10 <- paste("MCE.power", auc_vars, sep = ".")
names11 <- paste("av.time", auc_vars, sep = ".")
names12 <- paste("Missing_AUC", auc_vars, sep = ".")
names13 <- paste("cov_be", auc_vars, sep = ".")
names14 <- paste("MCE.cov_be", auc_vars, sep = ".")
colnames(average) <- c(names1, names6,names3, names8, names2, names7, names4, names9, names5, names10, names13, 
                       names14, names11, names12)
times <- grep("time", names(res), value = TRUE)
nsim <- max(res$sim)

for (i in 1:nrow(scenario)){
  for (j in 1:length(auc_vars)){
    
    meanAUC <- mean(res[which(res$scenario==i),auc_vars[j]], na.rm=TRUE)
    var <- sum((res[which(res$scenario==i), auc_vars[j]]-meanAUC)^2, na.rm=TRUE)
    average[i,j] <- mean(df[which(res$scenario==i),j], na.rm = TRUE) # mean bias 
    average[i,j+nvars] <- sqrt(var/(nsim*(nsim-1))) # MCE SE of bias
    average[i,j+nvars*2] <- sqrt(var/(nsim-1)) # empSE
    average[i,j+nvars*3] <- average[i,j+nvars*2]/sqrt(2*(nsim-1)) # MCE SE of empSE
    average[i,j+nvars*4] <- sqrt(mean(df[which(res$scenario==i),j+nvars], na.rm=TRUE)) # RMSE
    average[i,j+nvars*5] <- sqrt(sum(((df[which(res$scenario==i),j+nvars])-mean(df[which(res$scenario==i),j+nvars], na.rm = TRUE))^2, na.rm = TRUE)/(nsim*(nsim-1))) # MCE of MSE
    average[i,j+nvars*6] <- mean(df[which(res$scenario==i),j+nvars*2], na.rm=TRUE) # mean coverage
    average[i,j+nvars*7] <- sqrt(average[i,j+nvars*6]*(1-average[i,j+nvars*6])/nsim) # MCE SE of coverage
    average[i,j+nvars*8] <- mean(df[which(res$scenario==i),j+nvars*3], na.rm = TRUE) # power
    average[i,j+nvars*9] <- sqrt(average[i,j+nvars*8]*(1-average[i,j+nvars*8])/nsim) # MCE SE of power
    
    average[i,j+nvars*10] <- sum(1*(res[which(res$scenario==i),cil]<=meanAUC 
                                     & meanAUC<=res[which(res$scenario==i),ciu]), na.rm=TRUE)/nsim # mean bias-eliminated coverage
    average[i,j+nvars*11] <- sqrt(average[i,j+nvars*10]*(1-average[i,j+nvars*10])/nsim) # MCE SE of be-coverage
    
    if (auc_vars[j]=="AUC.IPL.LG") {
      average[i,j+nvars*12] <- mean(res[which(res$scenario==i),times[j]], na.rm=TRUE)/2
    } else if (auc_vars[j]=="AUC.IPL.NP") {
      average[i,j+nvars*12] <- mean(res[which(res$scenario==i),times[j-1]], na.rm=TRUE)/2
    } else {
      average[i,j+nvars*12] <- mean(res[which(res$scenario==i),times[j]], na.rm=TRUE) # time
    }
    average[i,j+nvars*13] <- sum(is.na(res[which(res$scenario==i),auc_vars[j]])) # missing
    
  }
}

scenario_res <- cbind(scenario,average)


