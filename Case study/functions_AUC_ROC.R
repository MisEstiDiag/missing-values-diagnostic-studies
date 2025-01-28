
##############      Functions for the simulations study     ####################
##############      Author: Katharina Stahlmann             ####################
##############      Date: 10.11.2023                        ####################


##############              Data generation                 ####################

# Description: function to generate simulated data
# parameters:
# N = sample size
# p = prevalence of target condition
# AUC_0 = true AUC
# korr = correlation between index test and covariates
# mech = missingness mechanism (MCAR, MAR, MNAR)
# pm = proportion of missing values
# mypattern = pattern of missing values (only in the index test)
# seed = seed for data generation

data_fun <- function(N, p, AUC_0, korr, mech, pm, mypattern, seed){
  sigma <- matrix(ncol=4, nrow=4, korr)
  diag(sigma) <- 1
  n1=N*p
  n0=N-n1
  mu=qnorm(AUC_0)*sqrt(2) # calculate mu of index test for diseased sample
  set.seed(seed)
  X0 <-as.data.frame(rmvnorm(n0, mean = c(0, 0, 5, 35), sigma = sigma,
                             method=c("eigen"), pre0.9_9994 = FALSE, checkSymmetry = TRUE)) # non-diseased population
  X1 <- as.data.frame(rmvnorm(n1, mean = c(mu, 0, 5, 35), sigma = sigma,
                              method=c("eigen"), pre0.9_9994 = FALSE, checkSymmetry = TRUE)) # diseased population; V1 is index test
  X0$D <- 0 # reference test non-diseased
  X1$D <- 1 # reference test diseased
  data <- rbind(X0, X1)
  data$V2 <- ifelse((data$V2<=0),1,2) # dichotomize V2
  # insert missing values with ampute function (mice package)
  amp <- ampute(data, prop = pm, patterns = mypattern, mech = mech)
  data_mis <- amp$amp
  return(data_mis)
}


##############                   Methods                    ####################

##############      MI2 by Long et al. 2011a                ####################


# Description: k-nearest neighbor multiple imputation based on a prediction and propensity score
# Output: AUC+CIs, Sensitivity values for false positive rate + CIs (ROC)

# parameter: 
# data = dataset (with missing values; full sample, columns are index test, 
#                 disease indicator, covariates)
# d = name of disease variable (binary) (must be 0=healthy, 1=diseased)
# ind = name of index test variable (continuous)
# cov_pred,cov_prop = list of covariates (as =c("C1","c2",..)); can be the same
#                    or different for the propensity (cov_prop) and prediction score (cov_pred)
# w = weight for prediction and propensity score model, can be altered
# kN = number of nearest neighbors, kL = number of imputation datasets
# tt = list of false positive rates
# required libraries: ks, pROC

##  main function (must be executed with the intern function PredPropScore)

#mi2 <- function(data, d, ind, cov_pred, cov_prop, w=c(0.5,0.5), kN=3, kL=10, tt, kLevel=0.95, seed) {
# ...
#}

# This method is proposed by Long et al. (2011a). Its code is included in the supplemental material of Cheng and Tang (2020). 


##############      MIB2 by Long et al. 2011a               ####################


# Description: k-nearest neighbor multiple imputation based on a prediction and propensity score plus bootstrap step
# Output: AUC+CIs, Sensitivity values for false positive rate + CIs (ROC)

# parameter: 
# data = dataset (with missing values; full sample, columns are index test, 
#                 disease indicator (must be 0=healthy, 1=diseased), covariates)
# d = name of disease variable (binary)
# ind = name of index test variable (continuous)
# cov_pred,cov_prop = list of covariates (as =c("C1","c2",..)); can be the same
#                       or different for the propensity (cov_prop) and prediction score (cov_pred)
# w = weight for prediction and propensity score model, can be altered
# kN = number of nearest neighbors, kL = number of imputation datasets
# tt = list of false positive rates
# required libraries: ks, pROC

##  main function (must be executed with the intern function PredPropScore (above))

#mi2boot <- function(data, d, ind, cov_pred, cov_prop, w=c(0.5,0.5), kN=3, kL=10, tt, kLevel, seed) {
# ...
#}

# This method is proposed by Long et al. (2011a). Its code is included in the supplemental material of Cheng and Tang (2020). 


##############      HDEL by Wang & Quin 2012;14                #################


# Description: Hot deck imputation with empirical likelihood
# Output: AUC+CIs, Sensitivity values for false positive rate + CIs (ROC)

# parameter: 
# data = dataset (with missing values; full sample, columns are index test, 
#                 disease indicator (must be 0=healthy, 1=diseased), covariates)
# d = name of disease variable (binary)
# ind = name of index test variable (continuous)
# a = alpha value (standard 0.05)
# tt = list of false positive rates
# l1, l2 = limits for confidence interval

## main function (must be executed with the two intern functions hatF and lambda)
hdel <- function(data, d, ind, a, l1, l2, tt, seed){
  #library(pROC)
  set.seed(seed)
  
  # data preparation
  col1 <- match(paste(ind),names(data)) # column number of index test
  col2 <- match(paste(d),names(data)) # column number of disease variable
  xx <- data[which(data[,col2]==1), col1] # only index test values of diseased sample
  yy <- data[which(data[,col2]==0), col1] # only index test values of healthy sample
  m <- length(xx)
  n <- length(yy)
  delta1 <- 1-is.na(xx)
  delta2 <- 1-is.na(yy)
  hatpi1=sum(delta1)/m; # missingness probability, here assumption that it is MCAR
  hatpi2=sum(delta2)/n;
  
  if (sum(!is.na(xx))==0 | sum(!is.na(yy))==0) {
    AUC.HDEL1=NA
    AUC.HDEL2=NA
    CIU.HDEL=NA
    CIL.HDEL=NA
    ROC.HDEL=NA
    ROC.CIL.HDEL=NA
    ROC.CIU.HDEL=NA
  } else {
    # hot deck imputation
    oxx <- xx[delta1==1] # only observed sample
    oyy <- yy[delta2==1]
    starX=sample(oxx,m,replace = TRUE, prob = NULL); # random hot deck imputation; draws random replacement for missing values from observed values
    starY=sample(oyy,n,replace = TRUE, prob = NULL); # random hot deck imputation
    impx=ifelse((is.na(xx)), starX, xx) # fill incomplete dataset with random hot deck draws
    impy=ifelse((is.na(yy)), starY, yy)
    
    ## calculate AUC
    AUC.HDEL1=mean(outer(impx,rep(1,length(impy)))>= outer(rep(1,length(impx)),impy) ); 
    # Wilcoxon Mann Whitney AUC, here >= instead of <= as in script, as we assume diseased persons have higher values than non-diseased persons
    data_HDEL <- as.data.frame(c(impx, impy))
    data_HDEL$D <- c(rep(1, m), rep(0, n))
    AUC.HDEL2=auc(data_HDEL$D, data_HDEL[,1])
    
    # calculate CI with empirical likelihood
    orderX=sort(impx);
    orderY=sort(impy);
    poolXY=c(orderX,orderY);
    Ri=(rank(poolXY))[1:m]; #rank function! not order
    Sj=(rank(poolXY))[(m+1):(m+n)];
    #estimate of s10
    s10=1/((m-1)*n^2) * ( sum((Ri-1:m)^2) - m * (mean(Ri)-(m+1)/2)^2 );
    #estimate of s01;
    s01=1/((n-1)*m^2) * ( sum((Sj-1:n)^2) - n * (mean(Sj)-(n+1)/2)^2 );
    #estimate of s
    s=(m*s01*(1-hatpi2+hatpi2^(-1)) + n*s10*(1-hatpi1+hatpi1^(-1)))/(m+n); 
    #EL
    U=1-hatF(impy,impx);
    w <- seq(l1, l2, 0.001) # how to select this region?, in values <0.8 the logr==NA...
    w <- sort(c(w, AUC.HDEL1)) # include calculated AUC in possible value list
    logr <- unlist(lapply(w, lam, U=U, s=0))
    #quantile
    qua=(m+n)*n*s / (m * sum( (1-U-AUC.HDEL1)^2 ) ) * qchisq(a,1); # quantile which defines the 95% confidence region
    ind <- which(logr<=qua)
    ci <- w[ind]
    #confidence interval
    CIU.HDEL=max(ci)
    CIL.HDEL=min(ci)
    
    ## ROC values (Sensitivities calculation)
    #library(ks)
    ROC.HDEL <- rep(0,length(tt))
    ROC.CIU.HDEL <- rep(0,length(tt))
    ROC.CIL.HDEL <- rep(0,length(tt))
    test <-  matrix(0, 81, length(tt))
    for (k in 1:length(tt)) {
      ROC.HDEL[k]=mean(impx >= quantile(impy,1-tt[k])); # estimated sensitivity given imputed data and given specificity
      if (ROC.HDEL[k]<1){
        #estimate of sigma
        sigma=ROC.HDEL[k]*(1-ROC.HDEL[k]);
        hatpi1=sum(delta1)/m;
        hatpi2=sum(delta2)/n;
        invFp=quantile(impy,1-tt[k]);
        den1=kde(x=impx, h=hpi(impx),eval.point=invFp); # kernel density estimates
        den2=kde(x=impy, h=hpi(impy),eval.point=invFp);
        denest1=den1$estimate;
        denest2=den2$estimate;
        sigma1=sigma*(1-hatpi1+1/hatpi1)+m/n*(1-tt[k])*(tt[k])*denest1^2/denest2^2*(1-hatpi2+1/hatpi2);
        # estimate CI
        w <- seq(0.2, 0.99, 0.01) 
        w <- sort(c(w, ROC.HDEL[k]))
        logr <- unlist(lapply(w, lam, U=U, s=tt[k]))
        test[,k] <- logr
        #quantile
        qua=sigma1/sigma*qchisq(a,1);
        ind <- which(logr<=qua & logr>0 )
        ci <- w[ind]
        #confidence interval
        ROC.CIU.HDEL[k]=max(ci)
        ROC.CIL.HDEL[k]=min(ci)
      }
      else{ROC.CIU.HDEL[k]=NA; ROC.CIL.HDEL[k]=NA}
    }
  }
  return(list(AUC.HDEL1=AUC.HDEL1, AUC.HDEL2=AUC.HDEL2, CIU.HDEL=CIU.HDEL, CIL.HDEL=CIL.HDEL, 
              ROC.HDEL=ROC.HDEL, ROC.CIU.HDEL=ROC.CIU.HDEL, ROC.CIL.HDEL=ROC.CIL.HDEL))
}

## intern functions
hatF<-function(Y,X)
{
  temp=1 * (outer(Y,rep(1,length(X)))<= outer(rep(1,length(Y)),X) ); # I think I heave to change the operator here too into <=
  Ef=apply(temp,MARGIN=1,FUN="mean");
  #return(temp)
  return(Ef);
}

lam <- function(w, U, s){
  if (s==0){W=1-U-w}
  else {W=1*(U<=s) - w}
  
  #solve lambda
  lambda=0;
  temp0<-1;
  for (j in 1:100)
  {
    der<- - mean(W^2/(1+temp0*W)^2);
    F<- mean(W/(1+temp0*W));
    temp1<-temp0-F/der;
    if (is.na(abs(temp1-temp0))){lambda<-temp0}
    else if (abs(temp1-temp0)<10^(-6)){
      lambda<-temp1;
      #print(lambda);
      break;
    }
    temp0<-temp1;
  }
  #log-EL ratio
  logr=2*sum(log(1+lambda*W));
  return(logr)
}



##############      AIPW by Long et al. 2011b                  #################


# Description: Augmented inverse probability weightin approach using a propsensity
            #  score model and a prediction score, plus bootstrap
# Output: AUC+CIs

# parameter: 
# data = dataset (with missing values; full sample, columns are index test, 
#                 disease indicator (must be 0=healthy, 1=diseased), covariates)
# d = name of disease variable (binary)
# ind = name of index test variable (continuous)
# cov_pred / cov_prop = list of covariates (as =c("C1","c2",..))

# main function (must be executed with intern functions prop & pred)

#aipw <- function(data, d, ind, cov_pred, cov_prop, seed){
# ...
#}

# This method is proposed by Long et al. (2011b). Its code is included in the supplemental material of Cheng and Tang (2020). 


##############      CONV by Bianco et al. 2023                 #################


# Description: Convolution-based estimators 
# Output: AUC, Sensitivity values for false positive rate (ROC)

# parameter: 
# data = dataset (with missing values; full sample, columns are index test, 
#                 disease indicator, covariates)
# d = name of disease variable (binary)
# ind = name of index test variable (continuous)
# cov_pred / cov_prop = list of covariates (as =c("C1","c2",..))
# tt = list of false positive rates

# main function (must be executed with the intern function roc_CONV)
#conv <- function(data, d, ind, cov_pred, cov_prop, tt){
 
 # ...

#}

# This method is described by Bianco et al. (2023). They provide the code to their method here: https://github.com/gboente/ROC-missing



##############      KER by Bianco et al. 2023                  #################


# Description: Kernel-based estimators 
# Output: AUC, Sensitivity values for false positive rate  (ROC)

# parameter: 
# data = dataset (with missing values; full sample, columns are index test, 
#                 disease indicator, covariates)
# d = name of disease variable (binary)
# ind = name of index test variable (continuous)
# cov_pred, cov_prop = list of covariates (as =c("C1","c2",..))
# tt = list of false positive rates

# main function (must be executed with the intern function roc_KER)
#ker <- function(data, d, ind, cov_pred, cov_prop, tt){
 # ...
#}
  
# This method is described by Bianco et al. (2023). They provide the code to their method here: https://github.com/gboente/ROC-missing



###########             IPL/NP by Cheng & Tang 2020             ################


# Description: parametric and non-parametric approach usind hybrid imputation 
            #  and empirical liklelihood
# Output: AUC+CIs, Sensitivity values for false positive rate + CIs (ROC)

# parameter: 
# data = dataset (with missing values; full sample, columns are index test, 
#                 disease indicator, covariates)
# d = name of disease variable (binary)
# ind = name of index test variable (continuous)
# cov = list of covariates (as =c("C1","c2",..))
# tt = list of false positive rates # tt=seq(0.01, 1-0.01, 0.01)
# kLevel = 1-alpha (set to 0.95 here)
# kK = number of imputations (set to 20 here)

#ipl <- function(data, d, ind, cov, tt, kLevel=0.95, kK=20, seed){
# ...
#}

# This method is described by Cheng and Tang (2020). Its code can be found in the supplemental material of their publication.


###########                   References                        ################

#Bianco AM, Boente G, González–Manteiga W, Pérez–González A. Estimators for ROC curves with missing biomarkers values and informative covariates. Statistical Methods & Applications. 2023.

#Cheng W, Tang N. Smoothed empirical likelihood inference for ROC curve in the presence of missing biomarker values. Biom J. 2020;62(4):1038-59.

#Long Q, Zhang X, Hsu C-H. Nonparametric multiple imputation for receiver operating characteristics analysis when some biomarker values are missing at random. Stat Med. 2011a;30(26):3149-61.

#Long Q, Zhang X, Johnson BA. Robust estimation of area under ROC curve using auxiliary variables in the presence of missing biomarker values. Biometrics. 2011b;67(2):559-67.

#Wang B, Qin G. Imputation-based empirical likelihood inference for the area under the ROC curve with missing data. Stat Interface. 2012;5(3):319-29.

#Wang B, Qin G. Empirical likelihood-based confidence intervals for the sensitivity of a continuous-scale diagnostic test with missing data. Commun Stat Theory Methods. 2014;43(15):3248-68.



