## ########################################################################## ##
###################     Simulation Study: Case study        ####################
##                          Author: KS                                        ##
##                      Date: November 2023                                   ##
##                        Out: AUC (ROC)                                      ##
## ########################################################################## ##


#####   Generate dataset with missing data    #####

amp <- ampute(data_full, prop = pm, patterns = mypattern, mech = mech)
summary(amp$amp)
data <- amp$amp

# Cave: for all functions, the reference test must take the values 0 and 1. The covariates should not be labelled factors, but numeric. 

# prepare amputed dataset
data <- data %>%
  mutate(
    ear = as.factor(ear),
    sitenum = as.numeric(sitenum),
    id = as_factor(id),
    gender = as.numeric(gender), #(1=F, 2=M)
    d = as.numeric(d),
    reftest = if_else(d==2, 1,
                      if_else(d==1, 0, NA))
  )


#####   calculate AUC and ROC by each method    #####

## CCA for true data
ci <- ci.auc(data_full$d,data_full$y1)
AUC.CCA_full <- ci[2]
CIU.CCA_full <- ci[3]
CIL.CCA_full <- ci[1]
# ROC
tt=seq(0.01, 1-0.01, 0.01)
roc <- roc(data_full$d, data_full$y1)
x <- rep(NA,length(tt))
for (i in 1:length(tt)){
  dist <- (tt[i]-(1-roc$specificities))^2
  x[i] <- order(dist)[1]
}
fpr <- roc$specificities[x]
ROC.CCA_full <- roc$sensitivities[x]
ci <- ci(roc, of = "se", specificities = 1-fpr)
ROC.CIL.CCA_full <- ci[,1]
ROC.CIU.CCA_full <- ci[,2]

## CCA for missing data
ci <- ci.auc(data$d,data$y1)
AUC.CCA <- ci[2]
CIU.CCA <- ci[3]
CIL.CCA <- ci[1]
# ROC
tt=seq(0.01, 1-0.01, 0.01)
roc <- roc(data$reftest, data$y1)
x <- rep(NA,length(tt))
for (i in 1:length(tt)){
  dist <- (tt[i]-(1-roc$specificities))^2
  x[i] <- order(dist)[1]
}
fpr <- roc$specificities[x]
ROC.CCA <- roc$sensitivities[x]
ci <- ci(roc, of = "se", specificities = 1-fpr)
ROC.CIL.CCA <- ci[,1]
ROC.CIU.CCA <- ci[,2]

## HDEL
res <- hdel(data=data, d="reftest", ind="y1", a=0.95, l1=0.5, l2=0.95, tt=seq(0.01, 1-0.01, 0.01), seed=548) # warnings indicate that <0.75 there is no logr
AUC.HDEL <- res$AUC.HDEL1
CIU.HDEL <- res$CIU.HDEL
CIL.HDEL <- res$CIL.HDEL
ROC.HDEL <- res$ROC.HDEL
ROC.CIL.HDEL <- res$ROC.CIL.HDEL
ROC.CIU.HDEL <- res$ROC.CIU.HDEL

## MI2
res <- mi2(data=data, d="reftest", ind="y1", cov_pred=c("currage","gender","sitenum"), cov_prop=c("currage","gender","sitenum"), 
           w=c(0.5,0.5), kN=3, kL=10, tt=seq(0.01, 1-0.01, 0.01), kLevel=0.95, seed=395)
AUC.MI2 <- res$AUC.MI2
CIU.MI2 <- res$CIU.MI2
CIL.MI2 <- res$CIL.MI2
ROC.MI2 <- res$ROC.MI2
ROC.CIL.MI2 <- res$ROC.CIL.MI2
ROC.CIU.MI2 <- res$ROC.CIU.MI2

## MIB2
res <- mi2boot(data=data, d="reftest", ind="y1", cov_pred=c("currage","gender","sitenum"), cov_prop=c("currage","gender","sitenum"), 
               w=c(0.5,0.5), kN=3, kL=10, tt=seq(0.01, 1-0.01, 0.01), kLevel=0.95, seed=482)
AUC.MIB2 <- res$AUC.MIB2
CIU.MIB2 <- res$CIU.MIB2
CIL.MIB2 <- res$CIL.MIB2
ROC.MIB2 <- res$ROC.MIB2
ROC.CIL.MIB2 <- res$ROC.CIL.MIB2
ROC.CIU.MIB2 <- res$ROC.CIU.MIB2

## mice
set.seed(954)
predMatrix <- make.predictorMatrix(data[,2:ncol(data)]) # in some cases predMatrix should be manually altered
impMethod <- make.method(data[,2:ncol(data)]) 
impMethod[] <- ""
idx <- names(which(colSums(is.na(data[,2:ncol(data)])) > 0))
impMethod[idx] <- "pmm" # pmm for all variables with missing data
imp <- mice(data[,2:ncol(data)], m=20, maxit=10, method = impMethod, predictorMatrix = predMatrix) # imputation
implist <- mitml::mids2mitml.list(imp)
# AUC
fit <- unlist(lapply(implist, function(d) auc(d$reftest, d$y1))) # estimate AUC
AUC.mice = mean(fit)
se <- unlist(lapply(implist, function(d) sqrt(var(auc(d$reftest, d$y1))))) # estimate CI
var.within = sum(se^2)/length(se) 
var.between = (sum((fit-AUC.mice)^2))/(length(fit)-1)
se.mice = sqrt(var.within+var.between+(var.between/length(se))) # pooled SE
CIU.mice = AUC.mice+1.96*se.mice # confidence interval
CIL.mice = AUC.mice-1.96*se.mice
# ROC
fit <- lapply(implist, function(d) roc(d$reftest, d$y1))
tt=seq(0.01, 1-0.01, 0.01)
# first ROC estimate for each imputed dataset and then pool
fprlist <- list()
roclist <- list()
for (j in 1:length(fit)) {
  
  x <- rep(NA,length(tt))
  for (i in 1:length(tt)){
    dist <- (tt[i]-(1-fit[[j]]$specificities))^2
    x[i] <- order(dist)[1]
  }
  fprlist[[j]] <- fit[[j]]$specificities[x]
  roclist[[j]] <- fit[[j]]$sensitivities[x]
  
}
fpr2 <- do.call(cbind.data.frame, fprlist)
fpr <- rowMeans(fpr2) # pooling across imputed datasets
roc <- do.call(cbind.data.frame, roclist)
ROC.mice <- rowMeans(roc)
ROC.CIL.mice <- rep(NA, 99)
ROC.CIU.mice <- rep(NA, 99)

## mix
data_mix <- data[,c(3,5,10,4,7:9)] # restructure, categorical variables first
data_mix$reftest[data_mix$reftest==0] <- 2 # 2=non-diseased, 1=diseased
s <- prelim.mix(data_mix, 3) # 1 gives the number of categorical variables in the dataset
thetahat <- em.mix(s) # get starting value
rngseed(1545) # set random number generator seed, must be set (otherwise da.mix wont work), change for simulation
implist=list()
for (m in 1:20) { # imputation
  theta <- da.mix(s,thetahat,steps=100) # data augmentation, thetahat provides the starting values
  implist[[m]] <- as.data.frame(imp.mix(s, theta, data_mix))
}
fit <- unlist(lapply(implist, function(d) auc(d$reftest, d$y1))) # estimate AUC
AUC.mix = mean(fit) 
se <- unlist(lapply(implist, function(d) sqrt(var(auc(d$reftest, d$y1))))) # estimate CI
var.within = sum(se^2)/length(se) # 20=no of imputation sets
var.between = (sum((fit-AUC.mix)^2))/(length(fit)-1)
se.mix = sqrt(var.within+var.between+(var.between/length(se))) # pooled SE
CIU.mix = AUC.mix+1.96*se.mix # confidence interval
CIL.mix = AUC.mix-1.96*se.mix
# ROC
fit <- lapply(implist, function(d) roc(d$reftest, d$y1))
tt=seq(0.01, 1-0.01, 0.01)
# first ROC estimate for each imputed dataset and then pool
fprlist <- list()
roclist <- list()
for (j in 1:length(fit)) {
  
  x <- rep(NA,length(tt))
  for (i in 1:length(tt)){
    dist <- (tt[i]-(1-fit[[j]]$specificities))^2
    x[i] <- order(dist)[1]
  }
  fprlist[[j]] <- fit[[j]]$specificities[x]
  roclist[[j]] <- fit[[j]]$sensitivities[x]
  
}
fpr2 <- do.call(cbind.data.frame, fprlist)
fpr <- rowMeans(fpr2) # pooling across imputed datasets
roc <- do.call(cbind.data.frame, roclist)
ROC.mix <- rowMeans(roc)
ROC.CIL.mix <- rep(NA, 99)
ROC.CIU.mix <- rep(NA, 99)

## mi
set.seed(275)
mdf <- missing_data.frame(data_mix)
mdf <- change(mdf, y = "y1", what = "imputation_method", to = "pmm")
imp <- mi(mdf, n.chains=20) # imputation
implist <- mi::complete(imp)
fit <- unlist(lapply(implist, function(d) auc(d$reftest, d$y1))) # estimate AUC
AUC.mi_ = mean(fit)
se <- unlist(lapply(implist, function(d) sqrt(var(auc(d$reftest, d$y1))))) # estimate CI
var.within = sum(se^2)/length(se) # 20=no of imputation sets
var.between = (sum((fit-AUC.mi_)^2))/(length(fit)-1)
se.mi = sqrt(var.within+var.between+(var.between/length(se))) # pooled SE
CIU.mi_ = AUC.mi_+1.96*se.mi # confidence interval
CIL.mi_ = AUC.mi_-1.96*se.mi
# ROC
fit <- lapply(implist, function(d) roc(d$reftest, d$y1))
tt=seq(0.01, 1-0.01, 0.01)
# first ROC estimate for each imputed dataset and then pool
fprlist <- list()
roclist <- list()
for (j in 1:length(fit)) {
  
  x <- rep(NA,length(tt))
  for (i in 1:length(tt)){
    dist <- (tt[i]-(1-fit[[j]]$specificities))^2
    x[i] <- order(dist)[1]
  }
  fprlist[[j]] <- fit[[j]]$specificities[x]
  roclist[[j]] <- fit[[j]]$sensitivities[x]
  
}
fpr2 <- do.call(cbind.data.frame, fprlist)
fpr <- rowMeans(fpr2) # pooling across imputed datasets
roc <- do.call(cbind.data.frame, roclist)
ROC.mi_ <- rowMeans(roc)
ROC.CIL.mi_ <- rep(NA, 99)
ROC.CIU.mi_ <- rep(NA, 99)

## Ker
# KER
res <- ker(data=data, d="reftest", ind="y1", cov_prop=c("currage","gender","sitenum"), cov_pred=c("currage","gender","sitenum"), tt=seq(0.01, 1-0.01, 0.01))
AUC.KER <- res$AUC.KER
CIU.KER <- res$CIU.KER
CIL.KER <- res$CIL.KER
ROC.KER <- res$ROC.KER
ROC.CIL.KER <- rep(NA, 99)
ROC.CIU.KER <- rep(NA, 99)

## AIPW
res <- aipw(data=data, d="reftest", ind="y1", cov_prop=c("currage","gender","sitenum"), cov_pred=c("currage","gender","sitenum"), seed=468)
AUC.AIPW <- res$AUC.AIPW
CIU.AIPW <- res$CIU.AIPW
CIL.AIPW <- res$CIL.AIPW

## conv (14 minutes)
t1 <- proc.time()
res <- conv(data=data, d="reftest", ind="y1", cov_prop=c("currage","gender","sitenum"), cov_pred=c("currage","gender","sitenum"), tt=seq(0.01, 1-0.01, 0.01))
t2 <- proc.time()
time <- t2-t1[3]
AUC.CONV <- res$AUC.CONV
CIU.CONV <- res$CIU.CONV
CIL.CONV <- res$CIL.CONV
ROC.CONV <- res$ROC.CONV
ROC.CIL.CONV <- rep(NA, 99)
ROC.CIU.CONV <- rep(NA, 99)



## dataset with results for AUC + CI
Method <- c("CCA full data", "CCA missing data", "HDEL", "MI2", "MIB2", "mice", "mix", "mi", "KER", "AIPW", "CONV")
AUC <- c(AUC.CCA_full, AUC.CCA, AUC.HDEL, AUC.MI2, AUC.MIB2, AUC.mice, AUC.mix, AUC.mi_, AUC.KER, AUC.AIPW, AUC.CONV)
CI_lower <- c(CIL.CCA_full, CIL.CCA, CIL.HDEL, CIL.MI2, CIL.MIB2, CIL.mice, CIL.mix, CIL.mi_, CIL.KER, CIL.AIPW, CIL.CONV)
CI_upper <- c(CIU.CCA_full, CIU.CCA, CIU.HDEL, CIU.MI2, CIU.MIB2, CIU.mice, CIU.mix, CIU.mi_, CIU.KER, CIU.AIPW, CIU.CONV)
res_auc <- as.data.frame(cbind(Method, AUC, CI_lower, CI_upper))



## dataset with ROC + CI
roc <- cbind(tt, ROC.CCA_full, ROC.CCA, ROC.HDEL, ROC.MI2, ROC.MIB2, ROC.mice, ROC.mix, ROC.mi_, ROC.KER, ROC.CONV)
roc_cil <- cbind(ROC.CIL.CCA_full, ROC.CIL.CCA, ROC.CIL.HDEL, ROC.CIL.MI2, ROC.CIL.MIB2, ROC.CIL.mice, ROC.CIL.mix, ROC.CIL.mi_, 
                 ROC.CIL.KER, ROC.CIL.CONV)
roc_ciu <- cbind(ROC.CIU.CCA_full, ROC.CIU.CCA, ROC.CIU.HDEL, ROC.CIU.MI2, ROC.CIU.MIB2, ROC.CIU.mice, ROC.CIU.mix, ROC.CIU.mi_,
                 ROC.CIU.KER, ROC.CIU.CONV)
roc <- as.data.frame(cbind(roc, roc_cil, roc_ciu))
