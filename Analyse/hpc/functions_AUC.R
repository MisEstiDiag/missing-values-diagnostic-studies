
##############      Functions for the simulations study     ####################
##############      Author: Katharina Stahlmann             ####################
##############      Date: 15.10.2023                        ####################


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
  mu=qnorm(AUC_0)*sqrt(2) # calculate mean of index test for diseased sample
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
# Output: AUC+CIs

# Parameter: 
# data = dataset (with missing values; full sample, columns are index test, 
#                 disease indicator, covariates)
# d = name of disease variable (binary) (must be 0=healthy, 1=diseased)
# ind = name of index test variable (continuous)
# cov_pred,cov_prop = list of covariates (as =c("C1","c2",..)); can be the same
#                     or different for the propensity (cov_prop) and prediction score variables (cov_pred)
# w = weight for prediction and propensity score model, can be altered
# kN = number of nearest neighbors 
# kL = number of imputation datasets
# tt = list of false positive rates
# required libraries: ks, pROC

##  main function (uses the intern function PredPropScore (see below))

mi2 <- function(data, d, ind, cov_pred, cov_prop, w=c(0.5,0.5), kN=3, kL=10, tt, kLevel=0.95, seed) {

  #library(pROC)
  set.seed(seed)
  
  # data preparation
  col1 <- match(paste(ind),names(data)) # column number of index test
  col2 <- match(paste(d),names(data)) # column number of disease variable
  data$delta <- 1-is.na(data[,col1]) # missingness indicator
  xx <- data[which(data[,col2]==1), ]
  yy <- data[which(data[,col2]==0), ]
  m <- dim(xx)[1]
  n <- dim(yy)[1]
  
  # if either xx or yy contain only missing values in the index test, then set AUC + CI to NA
  if (nrow(xx[which(xx$delta==1), ])==0 | nrow(yy[which(yy$delta==1), ])==0) {
    AUC.MI2=NA
    CIU.MI2=NA
    CIL.MI2=NA
  } else {
    # for diseased sample
    x.mi.pred <- PredPropScore(xx, kN=kN, MI=kL, cov_pred=cov_pred, cov_prop=cov_prop, col1=col1, w=w)
    # for non-diseased sample
    y.mi.pred <- PredPropScore(yy, kN=kN, MI=kL, cov_pred=cov_pred, cov_prop=cov_prop, col1=col1, w=w)
    
    ## AUC calculation
    # combine mi datasets
    data_MI2 <- as.data.frame(rbind(x.mi.pred, y.mi.pred))
    data_MI2$D <- c(rep(1, m), rep(0, n))
    
    AUC <- rep(0,kL) 
    se <- rep(0,kL)
    for (k in 1:kL){
      AUC[k]=auc(data_MI2$D, data_MI2[,k])
      se[k]=sqrt(var(auc(data_MI2$D, data_MI2[,k]))) #in pROC package: var(myROC) gives variance of AUC; square root of var(ROC) gives SD and SE (both!)
    }
    AUC.MI2=mean(AUC)
    var.within = sum(se^2)/kL
    var.between = (sum((AUC-AUC.MI2)^2))/(kL-1)
    se.MI2 = sqrt(var.within+var.between+(var.between/kL)) # pooled SE
    CIU.MI2=AUC.MI2+1.96*se.MI2 # confidence interval
    CIL.MI2=AUC.MI2-1.96*se.MI2
  }
  
  return(list(AUC.MI2=AUC.MI2, CIU.MI2=CIU.MI2, CIL.MI2=CIL.MI2))
}

##  intern function

# Parameters:
# xx = only diseased / non-diseased sample, respectively
# cov_pred,cov_prop = list of covariates (as =c("C1","c2",..)); can be the same
#                     or different for the propensity (cov_prop) and prediction score variables (cov_pred)
# w = weight for prediction and propensity score model, can be altered
# kN = number of nearest neighbors 
# MI = kL = number of imputation datasets
# col1 = column number of index test

PredPropScore <- function(xx, kN, MI, cov_pred, cov_prop, col1, w){
  # prediction score
  oxx <- xx[which(xx$delta==1), ] # only observed
  covars_pred <- paste(cov_pred, collapse = "+") # select covariates
  form <- paste("oxx[,col1] ~",covars_pred)
  alpha.x <- lm(as.formula(form), data=oxx)$coefficients
  df <- as.matrix(xx[ ,which((names(xx) %in% cov_pred)==TRUE)])
  x.pred <- cbind(1, df) %*% alpha.x
  x.pred.stand <- (x.pred - mean(x.pred)) / sd(x.pred) # standardize to mu=0, sd=1
    
  # propensity score
  covars_prop <- paste(cov_prop, collapse = "+") # select covariates
  form <- paste("xx$delta ~",covars_prop)
  beta.x <- glm(as.formula(form), binomial(link='logit'), 
                data=data.frame(xx))$coefficients
  df <- as.matrix(xx[ ,which((names(xx) %in% cov_prop)==TRUE)])
  x.prop <- cbind(1, df) %*% beta.x
  x.prop.stand <- (x.prop - mean(x.prop)) / sd(x.prop) # standardize to mu=0, sd=1
  
  # k-nearest neighbor selection
  loc.obs.x <- which(xx$delta == 1)
  m <- dim(xx)[1]
  mi.pred <- matrix(0, m, MI) # matrix with no of imputed datasets * no of observations of xx
  condition <- xx$delta==0
  for (k in (1:m)){ 
    if (condition[k]) { # for all missing observations of xx
      dist.x <- sqrt(w[1] * (x.pred.stand[k] - x.pred.stand)^2 + # distance to all observed obs for each missing obs separately
                         w[2] * (x.prop.stand[k] - x.prop.stand)^2)
      x.obs.loc <- intersect(order(dist.x), loc.obs.x)[1:kN] # list kn observed nearest neighbors
      loc.impu.x <- sample(x.obs.loc, MI, replace = TRUE, prob = NULL) # select one of the kn neighbors randomly for each MI imputed dataset
      mi.pred[k, ] <- xx[loc.impu.x, 1] # fill matrix lines of missing values with observed index test values values of the selected nearest neighbors
    } else { # for all observed values
      #sample(1:3, MI, replace = TRUE, prob = NULL) 
      mi.pred[k, ] <- xx[k, col1] # fill observed matrix values with observed values of index test 
    }
  }
  
  return(mi.pred=mi.pred)
}



##############      MIB2 by Long et al. 2011a               ####################


# Description: k-nearest neighbor multiple imputation based on a prediction and propensity score plus bootstrap step
# Output: AUC+CIs

# Parameter: 
# data = dataset (with missing values; full sample, columns are index test, 
#                 disease indicator, covariates)
# d = name of disease variable (binary) (must be 0=healthy, 1=diseased)
# ind = name of index test variable (continuous)
# cov_pred,cov_prop = list of covariates (as =c("C1","c2",..)); can be the same
#                     or different for the propensity (cov_prop) and prediction score (cov_pred)
# w = weight for prediction and propensity score model, can be altered
# kN = number of nearest neighbors
# kL = number of imputation datasets
# tt = list of false positive rates
# required libraries: ks, pROC

##  main function (uses the intern function PredPropScore (above))

mi2boot <- function(data, d, ind, cov_pred, cov_prop, w=c(0.5,0.5), kN=3, kL=10, tt, kLevel, seed) {
  #library(pROC)
  set.seed(seed)
  
  # data preparation
  col1 <- match(paste(ind),names(data)) # column number of index test
  col2 <- match(paste(d),names(data)) # column number of disease variable
  data$delta <- 1-is.na(data[,col1]) # missingness indicator
  xx <- data[which(data[,col2]==1), ]
  yy <- data[which(data[,col2]==0), ]
  
  # if either xx or yy contain only missing values in the index test, then set AUC + CI to NA
  if (nrow(xx[which(xx$delta==1), ])==0 | nrow(yy[which(yy$delta==1), ])==0) {
    AUC.MIB2=NA
    CIU.MIB2=NA
    CIL.MIB2=NA
  } else {
    # loop over mi datasets including bootstrap sample
    s <- length(tt)
    roc <- matrix(0, kL, s)
    var.roc <- matrix(0, kL, s)
    AUC <- rep(0,kL)
    se <- rep(0,kL)
    x.mi.pred <- NA
    y.mi.pred <- NA
    for (b in 1:kL){ 
      # for diseased sample
      m <- dim(xx)[1]
      xx.boot <- data.frame(matrix(ncol = ncol(xx), nrow = nrow(xx)))
      while (sum(is.na(xx.boot[,1]))==nrow(xx.boot)) {
        ind.x <- sample(1:m, m, replace = TRUE, prob = NULL) # bootstrap sample
        xx.boot <- xx[ind.x,]
      }
      x.mi.pred <- PredPropScore(xx=xx.boot, kN=kN, MI=1, cov_pred=cov_pred, cov_prop=cov_prop, col1=col1, w=w)
      
      # for the non-diseased sample
      n <- dim(yy)[1]
      yy.boot <- data.frame(matrix(ncol = ncol(yy), nrow = nrow(yy)))
      while (sum(is.na(yy.boot[,1]))==nrow(yy.boot)) {
        ind.y <- sample(1:n, n, replace = TRUE, prob = NULL) # bootstrap sample
        yy.boot <- yy[ind.y, ]
      }
      y.mi.pred <- PredPropScore(xx=yy.boot, kN=kN, MI=1, cov_pred=cov_pred, cov_prop=cov_prop, col1=col1, w=w)
      
      # AUC estimation
      data_MIB2 <- as.data.frame(c(x.mi.pred, y.mi.pred))
      data_MIB2$D <- c(rep(1, length(x.mi.pred)), rep(0, length(y.mi.pred)))
      AUC[b]=auc(data_MIB2$D, data_MIB2[,1])
      se[b]=sqrt(var(auc(data_MIB2$D, data_MIB2[,1]))) 
    }
    
    # Pool AUC
    AUC.MIB2 = mean(AUC)
    var.within = sum(se^2)/kL
    var.between = (sum((AUC-AUC.MIB2)^2))/(kL-1)
    se.MIB2 = sqrt(var.within+var.between+(var.between/kL)) # pooled SE
    CIU.MIB2=AUC.MIB2+1.96*se.MIB2 # confidence interval
    CIL.MIB2=AUC.MIB2-1.96*se.MIB2
  }
  
  return(list(AUC.MIB2=AUC.MIB2, CIU.MIB2=CIU.MIB2, CIL.MIB2=CIL.MIB2))
}



##############      HDEL by Wang & Quin 2012;14                #################


# Description: Hot deck imputation with empirical likelihood
# Output: AUC+CIs

# parameter: 
# data = dataset (with missing values; full sample, columns are index test, 
#                 disease indicator, covariates)
# d = name of disease variable (binary) (must be 0=healthy, 1=diseased)
# ind = name of index test variable (continuous)
# a = alpha value (standard 0.05)
# tt = list of false positive rates
# l1, l2 = limits for confidence interval

## main function (must be executed with the two intern functions hatF and lambda (see below))
hdel <- function(data, d, ind, a, l1, l2, tt, seed){
  # library(pROC)
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
    # Wilcoxon Mann Whitney AUC, here >= instead of <= as in original script, as we assume diseased persons have higher values than non-diseased persons
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
    w <- seq(l1, l2, 0.01) 
    w <- sort(c(w, AUC.HDEL1)) # include calculated AUC in possible value list
    logr <- unlist(lapply(w, lam, U=U, s=0))
    #quantile
    qua=(m+n)*n*s / (m * sum( (1-U-AUC.HDEL1)^2 ) ) * qchisq(a,1); # quantile which defines the 95% confidence region
    ind <- which(logr<=qua)
    ci <- w[ind]
    #confidence interval
    CIU.HDEL=max(ci)
    CIL.HDEL=min(ci)
  }
  
  return(list(AUC.HDEL1=AUC.HDEL1, AUC.HDEL2=AUC.HDEL2, CIU.HDEL=CIU.HDEL, CIL.HDEL=CIL.HDEL))
}

## intern functions
hatF<-function(Y,X)
{
  temp=1 * (outer(Y,rep(1,length(X)))<= outer(rep(1,length(Y)),X) )
  Ef=apply(temp,MARGIN=1,FUN="mean")
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


# Description: Augmented inverse probability weighting approach using a propsensity
#              score model and a prediction score, plus bootstrap
# Output: AUC+CIs

# parameter: 
# data = dataset (with missing values; full sample, columns are index test, 
#                 disease indicator, covariates)
# d = name of disease variable (binary) (must be 0=healthy, 1=diseased)
# ind = name of index test variable (continuous)
# cov_pred / cov_prop = list of covariates (as =c("C1","c2",..))

# main function (must be executed with intern functions prop & pred)

aipw <- function(data, d, ind, cov_pred, cov_prop, seed){
  
  set.seed(seed)
  
  # data preparation
  col1 <- match(paste(ind),names(data)) # column number of index test
  col2 <- match(paste(d),names(data)) # column number of disease variable
  x <- data[which(data[,col2]==0), col1] # only index test values of healthy sample
  y <- data[which(data[,col2]==1), col1] # only index test values of diseased sample
  zxprop <- data[which(data[,col2]==0), cov_prop]
  zyprop <- data[which(data[,col2]==1), cov_prop]
  deltax <- ifelse((is.na(x)), 0, 1) # missingness indicator
  deltay <- ifelse((is.na(y)), 0, 1)
  n1=length(x)
  n2=length(y)
  n=n1+n2
  yobs <- y[deltay==1] # only observed index test
  xobs <- x[deltax==1]
  
  zxpred <- data[which(data[,col2]==0), cov_pred] # only covariates
  zypred <- data[which(data[,col2]==1), cov_pred] 
  zxobs <- zxpred[which(deltax==1), cov_pred] # only covariates of subjects with observed index test values
  zyobs <- zypred[which(deltay==1), cov_pred] 
  ZX <-as.matrix(cbind(rep(1,n1), zxpred)) # matrix for calculating predicted values 
  ZY <-as.matrix(cbind(rep(1,n2), zypred))
  
  if (sum(!is.na(x))==0 | sum(!is.na(y))==0) {
    AUC.AIPW=NA
    CIU.AIPW=NA
    CIL.AIPW=NA
  } else {
    wbx <- prop(deltax, cov_prop, zxprop, n1)
    wby <- prop(deltay, cov_prop, zyprop, n2)
    predx <- pred(xobs, zxobs, cov_pred, ZX)
    predy <- pred(yobs, zyobs, cov_pred, ZY)
    sigbx <- predx$sig
    sigby <- predy$sig
    fitmbx <- predx$fitm
    fitmby <- predy$fitm
    
    XB<-x
    YB<-y
    deltaBX<-deltax
    deltaBY<-deltay
    
    # bootstrap step with replacement, stratified by disease status
    for(b in 1:200){ 
      xbobs <- NA
      while (sum(is.na(xbobs))==length(xbobs)) {
        indx<-sample(seq(1,n1), n1, replace = TRUE) # sample x values
        xb<-x[indx]
        deltabx<-deltax[indx]
        xbobs<-xb[deltabx==1] # only observed sample
      }
      zbxprop<-zxprop[indx, ]
      zbxpred<-zxpred[indx,]
      ZBX<-as.matrix(cbind(rep(1,n1), zbxpred))
      zbxobs<-zbxpred[which(deltabx==1),]
      XB<-cbind(XB,xb)
      deltaBX<-cbind(deltaBX,deltabx)
      
      ybobs <- NA
      while (sum(is.na(ybobs))==length(ybobs)) {
        indy<-sample(seq(1,n2),n2,replace = TRUE) # sample obs numbers of y values
        yb<-y[indy] # select y, r and z values for sampled obs numbers
        deltaby<-deltay[indy]
        ybobs<-yb[deltaby==1] # only observed sample
      }
      
      zbyprop<-zyprop[indy, ]
      zbypred<-zypred[indy,]
      ZBY<-as.matrix(cbind(rep(1,n2), zbypred))
      zbyobs<-zbypred[which(deltaby==1),]
      YB<-cbind(YB,yb)
      deltaBY<-cbind(deltaBY,deltaby)
      
      # IPW/pred for x
      wwx <- prop(deltabx, cov_prop, zbxprop, n1)
      wbx <- cbind(wbx, wwx)
      predx <- pred(xbobs, zbxobs, cov_pred, ZBX)
      fitmbx <- cbind(fitmbx, predx$fitm)
      sigbx <- c(sigbx, predx$sig)
      
      # IPW/pred for y
      wwy <- prop(deltaby, cov_prop, zbyprop, n2)
      wby <- cbind(wby, wwy)
      predy <- pred(ybobs, zbyobs, cov_pred, ZBY)
      fitmby <- cbind(fitmby, predy$fitm)
      sigby <- c(sigby, predy$sig)
    }
    
    # calculate theta across bootstrap samples
    thetab <- 0
    theta <- rep(0, length(sigbx))
    NN <- rep(0, length(sigbx))
    xx <- ifelse((is.na(XB)),0,XB)
    yy <- ifelse((is.na(YB)),0,YB)
    for (i in 1:n1){
      for (j in 1:n2){
        theta<-theta+deltaBY[j,]*deltaBX[i,]*(yy[j,]>xx[i,])/n1/n2/wby[j,]/wbx[i,]-(deltaBY[j,]*deltaBX[i,]-wby[j,]*wbx[i,])/(wby[j,]*wbx[i,])*pnorm((fitmby[j,]-fitmbx[i,])/(sqrt(sigbx^2+sigby^2)),0,1)/n1/n2
        NN<-NN+deltaBY[j,]*deltaBX[i,]/n1/n2/wby[j,]/wbx[i,]
      }
    }
    thetab <- theta/NN
    AUC.AIPW <- thetab[1]
    CIU.AIPW <- AUC.AIPW+1.96*sqrt(var(thetab, na.rm=T))
    CIL.AIPW <- AUC.AIPW-1.96*sqrt(var(thetab, na.rm=T))
  }
  
  return(list(AUC.AIPW=AUC.AIPW, CIU.AIPW=CIU.AIPW, CIL.AIPW=CIL.AIPW))
}

# intern functions

# Parameters:
# delta = missingness indicator of only diseased (or non-diseased) sample
# cov = covariates for propensity score
# z = dataset containing the covariates with only diseased (or non-diseased) subjects
# n = sample size of diseased (or non-diseased) sample

prop <- function(delta, cov, z, n){
  covars <- paste(cov, collapse = "+") # select covariates
  form <- paste("delta ~",covars)
  logfit <- glm(as.formula(form), data=z, family = "binomial")
  w<-predict(logfit,type="response")
  w<-w*sum(delta/w)/n
  return(w)
}

# Parameters:
# obs = only observed index test values of diseased (or non-diseased) sample
# zobs = dataset with covariates of subjects with observed index test values of diseased (or non-diseased) sample
# cov = covariates for prediction score
# Z = matrix for calculating predicted values (including covariates for all diseased (or non-diseased) subjects)

pred <- function(obs, zobs, cov, Z){
  covars <- paste(cov, collapse = "+") # select covariates
  form <- paste("obs ~",covars)
    lmfit <- lm(as.formula(form), data=zobs)
    coef <- lmfit$coefficients
    fitm <- Z%*%coef
    sig <- summary(lmfit)$sigma
    return(list(fitm=fitm, sig=sig))
}



##############      CONV by Bianco et al. (2023)               #################


# Description: Convolution-based estimators 
# Output: AUC, Sensitivity values for false positive rate (ROC)

# Parameter: 
# data = dataset (with missing values; full sample, columns are index test, 
#                 disease indicator, covariates)
# d = name of disease variable (binary)
# ind = name of index test variable (continuous)
# cov_pred, cov_prop = list of covariates (as =c("C1","c2",..))
# tt = list of false positive rates

# main function (must be executed with the intern function roc_CONV)

#conv <- function(data, d, ind, cov_pred, cov_prop, tt){
 #...
#}

# This method is described by Bianco et al. (2023). They provide the code to their method here: https://github.com/gboente/ROC-missing



##############      KER by Bianco et al. (2023)                #################


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
#  ...
#}

# This method is described by Bianco et al. (2023). They provide the code to their method here: https://github.com/gboente/ROC-missing



###########             IPL/NP by Cheng & Tang 2020             ################


# Description: parametric and non-parametric approach usind hybrid imputation 
#              and empirical liklelihood
# Output: AUC+CIs, Sensitivity values for false positive rate + CIs (ROC)

# Parameter: 
# data = dataset (with missing values; full sample, columns are index test, 
#                 disease indicator, covariates)
# d = name of disease variable (binary)
# ind = name of index test variable (continuous)
# cov = list of covariates (as =c("C1","c2",..))
# tt = list of false positive rates
# kLevel = 1-alpha (set to 0.95 here)
# kK = number of imputations (set to 20 here)

# main function (employs several auxiliary functions (see below))

ipl <- function(data, d, ind, cov, tt, kLevel=0.95, kK=20, seed){
  # packages
  #library(mvtnorm);
  #library(KernSmooth);
  #library(rgl); library(ks);
  set.seed(seed) 
  
  # data preparation
  col1 <- match(paste(ind),names(data)) # column number of index test
  col2 <- match(paste(d),names(data)) # column number of disease variable
  delta <- 1-is.na(data[,col1])
  z <- as.data.frame(data[ ,cov]) 
  coll <- 1:ncol(z)
  AIC <- rep(0,length(coll))
  for( i in coll){ # select covariate
    AIC[i] <- glm(delta ~ z[,i], binomial(link='logit'))$aic	# association between possible covariate and missingness indicator, extract aic 
  }
  ii <- which.min(AIC)  # select covariate with lowest AIC, + 1 to assign the right column in dataset
  dat <- cbind(delta, z[,ii], data[,c(col1,col2)])
  if (length(which(is.na(dat[, 2]))) != 0) dat <- dat[-which(is.na(dat[, 2])), ] # delete obs with missing values in the covariate
  dat[which(is.na(dat[, 3])), 3] <- 0  # set all NAs to 0 in index test
  dat[, 2] <- (dat[, 2] / max(dat[, 2])) * 10 
  dat[, 3] <- (dat[, 3] / max(dat[, 3])) * 10 
  xx <- dat[which(dat[,4]==1), c(1:3)] # the diseased datasets
  yy <- dat[which(dat[,4]==0), c(1:3)] # the non-diseased datasets
  
  # paramerter
  critvalue <- qchisq(kLevel, 1)
  s <- length(tt)
  crit.ka <- c(4.956850, 4.828210, 4.849716, 4.996722)
  crit.mi <- c(4.890063, 4.962257, 5.135554, 5.033429)
  c1=1; c2=1; c3=1
  
  # Analysis
  m <- dim(xx)[1]
  n <- dim(yy)[1]  
  s <- length(tt)
  roc.ip.logit <- rep(0, s)
  low.ip.logit <- rep(0, s)
  upp.ip.logit <- rep(0, s)
  Lvalue.ip.logit <- rep(0, s)
  Uvalue.ip.logit <- rep(0, s)
  roc.ip.nopa <- rep(0, s)
  low.ip.nopa <- rep(0, s)
  upp.ip.nopa <- rep(0, s)
  Lvalue.ip.nopa <- rep(0, s)
  Uvalue.ip.nopa <- rep(0, s)
  
  ########## 'IP' method using logistic missing mechanism #########
  res.pz.logit <- MissingRateLogit(xx, yy)
  pz1.logit <- res.pz.logit$pz1
  pz2.logit <- res.pz.logit$pz2
  
  # MI part
  # imputation value matrix of the missing biomarker by MI method 
  z1 <- xx[, 2]
  z2 <- yy[, 2]
  h.z1 <- c2 * sd(z1) * (m^(- 1 / 3)) # bandwidth
  h.z2 <- c2 * sd(z2) * (n^(- 1 / 3))
  res.mi <- ImpuMIIP(xx, yy, kK, h.z1, h.z2) # imputed datasets
  
  x.mi <- res.mi$x.mi # matrix mit imputed index test values 
  y.mi <- res.mi$y.mi
  
  ox <- xx[xx[, 1]==1, 3]
  oy <- yy[yy[, 1]==1, 3]
  h.ox <- c1 * sd(ox) * (length(ox)^(- 1 / 3)) # bandwidth
  h.oy <- c1 * sd(oy) * (length(oy)^(- 1 / 3))
  
  for (i in 1:s){
    cat('IP-logit', i, "\n")
    roc.ip.logit[i] <- optimize(function(theta) ElrForThetaIP(theta, 
                                                              xx, yy, tt[i], x.mi, y.mi, pz1.logit, pz2.logit, h.ox, h.oy)$"-2LLR",
                                c(0, 1), tol = 1e-8)$minimum      
    ul.roc <- findUL(step=0.2, initStep=0, fun=ElrForThetaIP, 
                     MLE=roc.ip.logit[i], level=critvalue, xx=xx, yy=yy, 
                     para=tt[i], x.mi=x.mi, y.mi=y.mi, pz1=pz1.logit, pz2=pz2.logit,
                     h.ox=h.ox, h.oy=h.oy)
    low.ip.logit[i] <- ul.roc$Low
    upp.ip.logit[i] <- ul.roc$Up
    Lvalue.ip.logit[i] <- ul.roc$Lvalue
    Uvalue.ip.logit[i] <- ul.roc$Uvalue
  }
  
  # 'IP' method using nonparametric missing mechanism 
  h.z1.nopa <- c3 * sd(z1) * (m^(- 1 / 3))
  h.z2.nopa <- c3 * sd(z2) * (n^(- 1 / 3))
  res.pz.nopa <- MissingRateNopa(xx, yy, h.z1.nopa, h.z2.nopa)
  pz1.nopa <- res.pz.nopa$pz1
  pz2.nopa <- res.pz.nopa$pz2
  for (i in 1:s){ 
    cat('IP-nopa', i, "\n")
    roc.ip.nopa[i] <- optimize(function(theta) ElrForThetaIP(theta, xx,
                                                             yy, tt[i], x.mi, y.mi, pz1.nopa, pz2.nopa, h.ox, h.oy)$"-2LLR", 
                               c(0, 1), tol = 1e-8)$minimum      
    ul.roc <- findUL(step=0.2, initStep=0, fun=ElrForThetaIP, 
                     MLE=roc.ip.nopa[i], level=critvalue, xx=xx, yy=yy, 
                     para=tt[i], x.mi=x.mi, y.mi=y.mi, pz1=pz1.nopa, pz2=pz2.nopa,
                     h.ox=h.ox, h.oy=h.oy)
    low.ip.nopa[i] <- ul.roc$Low
    upp.ip.nopa[i] <- ul.roc$Up
    Lvalue.ip.nopa[i] <- ul.roc$Lvalue
    Uvalue.ip.nopa[i] <- ul.roc$Uvalue
  }  
  
  # calculate AUC
  AUC.IPL.LG=mean(roc.ip.logit)
  AUC.IPL.NP=mean(roc.ip.nopa)
  
  # calculate CI LG
  q0=AUC.IPL.LG*(1-AUC.IPL.LG)
  q1=(AUC.IPL.LG/(2-AUC.IPL.LG))-AUC.IPL.LG^2
  q2=(2*AUC.IPL.LG^2)/(1+AUC.IPL.LG)-AUC.IPL.LG^2
  CIU.IPL.LG=AUC.IPL.LG+1.96*sqrt((q0+(m-1)*q1+(n-1)*q2)/(n*m))
  CIL.IPL.LG=AUC.IPL.LG-1.96*sqrt((q0+(m-1)*q1+(n-1)*q2)/(n*m))
  
  # calculate CI NP
  q0=AUC.IPL.NP*(1-AUC.IPL.NP)
  q1=(AUC.IPL.NP/(2-AUC.IPL.NP))-AUC.IPL.NP^2
  q2=(2*AUC.IPL.NP^2)/(1+AUC.IPL.NP)-AUC.IPL.NP^2
  CIU.IPL.NP=AUC.IPL.NP+1.96*sqrt((q0+(m-1)*q1+(n-1)*q2)/(n*m))
  CIL.IPL.NP=AUC.IPL.NP-1.96*sqrt((q0+(m-1)*q1+(n-1)*q2)/(n*m))
  
  res <- list()                
  res$roc <- rbind(roc.ip.logit, roc.ip.nopa)
  res$low <- rbind(low.ip.logit, low.ip.nopa)
  res$upp <- rbind(upp.ip.logit, upp.ip.nopa)
  res$auc <- c(AUC.IPL.LG, AUC.IPL.NP)
  res$auclow <- c(CIL.IPL.LG, CIL.IPL.NP)
  res$aucupp <- c(CIU.IPL.LG, CIU.IPL.NP)
  
  return(list(res=res, 
              L.ip.logit=Lvalue.ip.logit, U.ip.logit=Uvalue.ip.logit,
              L.ip.nopa=Lvalue.ip.nopa, U.ip.nopa=Uvalue.ip.nopa))
}


# intern functions

# Let missing mechanism be the logistic model
MissingRateLogit <- function(xx, yy){
  # estimate response rate for diseased population
  xx <- data.frame(xx)      
  lgst.x <- glm(xx[, 1]~xx[, 2], binomial(link='logit'), data=xx)
  pz1 <- predict(lgst.x, type='response') 
  pz1[which(pz1 < 0.05)] <- 0.05
  # estimate response rate for diseased population
  yy <- data.frame(yy)      
  lgst.y<-glm(yy[, 1]~yy[, 2], binomial(link='logit'), data=yy)
  pz2 <- predict(lgst.y, type='response') 
  pz2[which(pz2 < 0.05)] <- 0.05
  a <- ls()
  rm(list=a[which(a != 'pz1' & a != 'pz2')]) 
  list("pz1" = pz1, "pz2" = pz2)
}

# imputation value matrix of the missing biomarker by MI method 
ImpuMIIP <- function(xx, yy, kK, h.z1, h.z2){
  m <- dim(xx)[1]
  ox <- xx[xx[, 1] == 1, 3]
  x.mi <- matrix(0, m, kK)
  for (k in 1:m) {
    px <- FzMI(xx[k, 2], xx, h.z1)
    # the imputated kK biomarkers of the k-th biomarker are stored the k-th row  
    # of imputated matrix "x.mi".
    x.mi[k, ] <- sample(ox, kK, replace = T, prob = px)         
  }   
  n <- dim(yy)[1]
  oy <- yy[yy[, 1] == 1, 3]
  y.mi <- matrix(0, n, kK)
  for (k in 1:n){
    py <- FzMI(yy[k, 2], yy, h.z2)
    y.mi[k, ] <- sample(oy, kK, replace = T, prob = py)       
  } 
  a <- ls()
  rm(list=a[which(a != 'x.mi' & a != 'y.mi')]) 
  list("x.mi"=x.mi, "y.mi"=y.mi)
}

# a kernel estimate of p(x|zi)
FzMI <- function(zz, xx, h.z){
  TINY <- sqrt( .Machine$double.xmin )
  delta <- xx[, 1]
  oz <- xx[delta == 1, 2]
  weight <- dnorm( (oz - zz) / h.z)
  if (max(weight) < TINY) weight <- weight + TINY
  weight.all <- sum(weight)
  p <- weight / weight.all  
  a <- ls()
  rm(list=a[which(a != 'p')])
  return(p) 
}

# moment function by hybrid imputation (IP) method 
MomentIP <- function(xx, yy, quan, x.mi, y.mi, pz1, pz2, h.ox, h.oy){
  m <- dim(xx)[1]
  delta1 <- xx[, 1]
  x <- xx[, 3]   
  n <- dim(yy)[1]
  delta2 <- yy[, 1]
  y <- yy[, 3]
  # smoothed moment   
  mm1 <- pnorm((quan - x) / h.ox)
  mm2 <- pnorm((quan - y) / h.oy)
  mm1_mi <- rowMeans(pnorm((quan - x.mi) / h.ox))  
  mm2_mi <- rowMeans(pnorm((quan - y.mi) / h.oy))	   
  mm1.f <- delta1 / pz1 * mm1 + (1 - delta1 / pz1) * mm1_mi 
  mm2.f <- delta2 / pz2 * mm2 + (1 - delta2 / pz2) * mm2_mi  
  a <- ls()
  rm(list=a[which(a != 'mm1.f' & a != 'mm2.f')]) 
  list("mm1.f" = mm1.f, "mm2.f" = mm2.f)
}

# Find lambda given moment in terms of sum(m/(1+lam*m))=0
SeekingLambda <- function( x, maxit=25, gradtol=1e-7, 
                           svdtol = 1e-9, itertrace=FALSE ){
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  if( n <= p )
    stop("Need more observations.")
  
  z <- t( t(x) )
  #
  #  Scale the problem, by a measure of the size of a 
  # typical observation.  Add a tiny quantity to protect
  # against dividing by zero in scaling.  Since z*lam is
  # dimensionless, lam must be scaled inversely to z.
  #
  TINY <- sqrt( .Machine$double.xmin )
  scale <- mean( abs(z) ) + TINY
  z <- z/scale
  lam <- rep(0,p)
  #
  #     Take some precaution against users specifying
  # tolerances too small.
  #
  if( svdtol < TINY )svdtol <- TINY
  if( gradtol < TINY)gradtol <- TINY
  #
  #    Preset the weights for combining Newton and gradient
  # steps at each of 16 inner iterations, starting with
  # the Newton step and progressing towards shorter vectors
  # in the gradient direction.  Most commonly only the Newton
  # step is actually taken, though occasional step reductions
  # do occur.
  #
  nwts <- c( 3^-c(0:3), rep(0,12) )
  gwts <- 2^( -c(0:(length(nwts)-1)))
  gwts <- (gwts^2 - nwts^2)^.5
  gwts[12:16] <- gwts[12:16] * 10^-c(1:5)
  #
  #    Iterate, finding the Newton and gradient steps, and
  # choosing a step that reduces the objective if possible.
  #
  nits <- 0
  gsize <- gradtol + 1
  while(  nits<maxit && gsize > gradtol ){
    arg  <- 1 + z %*% lam
    wts1 <- as.vector( llogp(arg, 1/n) )
    wts2 <- as.vector( -llogpp(arg, 1/n) )^.5
    grad <- as.matrix( -z*wts1 )
    #############grad <- as.vector( apply( grad, 2, sum ) )
    grad <- as.vector(rowsum(grad, rep(1, nrow(grad)) ) )
    gsize <- mean( abs(grad) )
    hess <- z*wts2
    #                                   -1
    #    The Newton step is -(hess'hess)    grad,
    #  where the matrix hess is a sqrt of the Hessian.
    #  Use svd on hess to get a stable solution.
    #
    
    ## may try La.svd() in R (v. > 1.0) for better LAPACK.
    ## or use QR decomposition on hess to solve it.
    svdh <- svd( hess )
    ##  svdh <- La.svd( hess )
    if( min(svdh$d) < max(svdh$d)*svdtol )
      svdh$d <- svdh$d + max(svdh$d)*svdtol
    nstep <- svdh$v %*% (t(svdh$u)/svdh$d)
    ## nstep <- t(svdh$vt) %*% (t(svdh$u)/svdh$d)
    nstep <- as.vector( nstep %*% matrix(wts1/wts2,n,1) )
    gstep <- -grad
    if(  sum(nstep^2) < sum(gstep^2) )
      gstep <- gstep*(sum(nstep^2)^.5/sum(gstep^2)^.5)
    ologelr <- -sum( llog(arg,1/n) )
    ninner <- 0
    for(  i in 1:length(nwts) ){
      nlogelr <- logelr( z, lam+nwts[i]*nstep+gwts[i]*gstep )
      if( nlogelr < ologelr ){
        lam <- lam+nwts[i]*nstep+gwts[i]*gstep
        ninner <- i
        break
      }
    }
    nits <- nits+1
    if(  ninner==0  )nits <- maxit
    # if( itertrace )
    # cat('# ',nits,': lambda = (',lam,'), logelr = (',nlogelr,'), 
    # gsize = (',gsize,'), ninner = ',ninner, '\n', sep = '')
  }
  
  list(lambda = lam/scale, "-2LLR" = -2*nlogelr, Pval = 1-pchisq(-2*nlogelr, 
      df=p), grad=grad*scale, hess=t(hess)%*%hess*scale^2, wts=wts1, nits=nits)
}

#  A variation with a second knot at a large value
#  M did not appear to work as well.
#
# Note1: The function llog() is equal to the natural logarithm on the interval 
#   from eps >0 to infinity.Between -infinity and eps, llog() is a quadratic.
# Note2: Llog() are continuous across the "knot" at eps.
# Note3: The cutoff point, eps, is usually 1/n, where n is the number of 
#   observations.  Unless n is extraordinarily large, dividing by eps is not 
#   expected to cause numerical difficulty.
# Note4: A variation with a second knot at a large value M did not appear to 
#   work as well.

llog <- function(z, eps ){
  ans <- z
  avoidNA <- !is.na(z)
  lo <- (z<eps) & avoidNA  
  ans[ lo  ] <- log(eps) - 1.5 + 2*z[lo]/eps - 0.5*(z[lo]/eps)^2
  ans[ !lo ] <- log( z[!lo] )
  ans
}

# Note1:llogp() and llogpp() are the first derivatives of Llog(). 
# Note2:llogp() are continuous across the "knot" at eps.
llogp <- function( z, eps ){
  ans <- z
  avoidNA <- !is.na(z)    ###added 3/2012
  lo <- (z<eps) & avoidNA
  ans[ lo  ] <- 2.0/eps - z[lo]/eps^2
  ans[ !lo ] <- 1/z[!lo]
  ans
}

# Note1:llogpp() and llogpp() are the two derivatives of Llog(). 
# Note2:llogpp() are continuous across the "knot" at eps.
llogpp <- function( z, eps ){
  ans <- z
  avoidNA <- !is.na(z) 
  lo <- (z<eps) & avoidNA    ### added same avoidNA as above
  ans[ lo  ] <- -1.0/eps^2
  ans[ !lo ] <- -1.0/z[!lo]^2
  ans
}

# The log of el ratio statistic 
logelr <- function(g, lam ){ 
  g <- as.matrix(g)
  n <- nrow(g)
  p <- ncol(g)
  if(  n <= p  )
    stop("Need more observations than variables in logelr.")
  arg <- 1 + g %*% lam
  return( - sum( llog(arg,1/n) ) ) 
}

# optimize quantile  fixed theta by Hybrid imputation (IP)
ElrForQuanIP <- function(quan, xx, yy, para, x.mi, y.mi, pz1, pz2, h.ox, h.oy, maxit=25, 
                         gradtol=1e-7, svdtol = 1e-9, itertrace=FALSE){
  theta <- para[1]
  ttt <- para[2]  
  res.mm <- MomentIP(xx, yy, quan, x.mi, y.mi, pz1, pz2, h.ox, h.oy)
  moment1 <- res.mm$mm1.f - (1 - theta) 
  moment2 <- res.mm$mm2.f - (1 - ttt) 
  lambda1 <- SeekingLambda(moment1,maxit = maxit, gradtol = gradtol, 
                           svdtol = svdtol, itertrace=itertrace)$lambda
  lambda2 <- SeekingLambda(moment2,maxit = maxit, gradtol = gradtol, 
                           svdtol = svdtol, itertrace=itertrace)$lambda
  elr <- -2 * (logelr(moment1, lambda1) + logelr(moment2, lambda2))
  if( itertrace )
    cat('# quan:', quan, ': -2LLR = (', elr, '), lambda1 = (', lambda1, '), 
      lambda2 = ', lambda2, '\n')
  a <- ls()
  rm(list=a[which(a != 'elr' & a != 'lambda1' & a !='lambda2')])
  list("-2LLR"=elr,"lambda1" = lambda1, "lambda2" = lambda2)
}


# optimize theta by Hybrid imputation (IP)
ElrForThetaIP <- function(theta, xx, yy, para, x.mi, y.mi, pz1, pz2, h.ox, h.oy, maxit=25, 
                          gradtol=1e-7, svdtol = 1e-9, itertrace=FALSE){
  oy <- yy[yy[, 1]==1, 3]
  a <- min(oy)
  b <- max(oy)
  ttt <- para 
  pa <- c(theta, ttt)                        
  quan <- optimize(function(quan) ElrForQuanIP(quan, xx, yy, pa, x.mi, y.mi, 
                                               pz1, pz2, h.ox, h.oy, maxit = maxit, gradtol = gradtol, svdtol = svdtol, 
                                               itertrace=itertrace)$"-2LLR", c(a-3, b+3), tol = 1e-8)$minimum 
  res.mm <- MomentIP(xx, yy, quan, x.mi, y.mi, pz1, pz2, h.ox, h.oy)                     
  moment1 <-  res.mm$mm1.f- (1 - theta) 
  moment2 <-  res.mm$mm2.f- (1 - ttt) 
  lambda1 <- SeekingLambda(moment1,maxit = maxit, gradtol = gradtol, 
                           svdtol = svdtol, itertrace=itertrace)$lambda
  lambda2 <- SeekingLambda(moment2,maxit = maxit, gradtol = gradtol, 
                           svdtol = svdtol, itertrace=itertrace)$lambda
  elr <- -2 * (logelr(moment1, lambda1) + logelr(moment2, lambda2))
  if( itertrace )
    cat('# theta:',theta,'quan:',quan,':,-2LLR=(',elr,'),lambda1=(',lambda1,'),
     lambda2=',lambda2,'\n')    
  a <- ls()
  rm(list=a[which(a != 'elr' & a != 'theta' & a !='quan' & a != 'lambda1' 
                  & a !='lambda2')])
  list("-2LLR"=elr, "theta"=theta, "quan"=quan,"lambda1"=lambda1, 
       "lambda2"=lambda2)
}

# Find lower and upper bound
findUL <- function (step = 0.01, initStep = 0, fun, MLE, level = 3.84, 
                    ...) {
  value <- 0
  step1 <- step
  Lbeta <- MLE - initStep
  for (i in 1:8) {
    while (value < level) {
      Lbeta <- Lbeta - step1
      value <- fun(Lbeta, ...)$"-2LLR"
    }
    Lbeta <- Lbeta + step1
    step1 <- step1/10
    value <- fun(Lbeta, ...)$"-2LLR"
  }
  value1 <- value
  value <- 0
  Ubeta <- MLE + initStep
  for (i in 1:8) {
    while (value < level) {
      Ubeta <- Ubeta + step
      value <- fun(Ubeta, ...)$"-2LLR"
    }
    Ubeta <- Ubeta - step
    step <- step/10
    value <- fun(Ubeta, ...)$"-2LLR"
  }
  if ((value1>level) | (value>level)) warning("Something wrong. Check the MLE 
      and step inputs.")
  a <- ls()
  rm(list=a[which(a != 'Lbeta' & a != 'Ubeta' & a !='step1' & a != 'step' 
                  & a !='value1' & a !='value')]) 
  return(list(Low = Lbeta, Up = Ubeta, FstepL = step1, FstepU = step, 
              Lvalue = value1, Uvalue = value))
}


# Let missing mechanism be nonparametric model
MissingRateNopa <- function(xx, yy, h.z1, h.z2){
  # estimate response rate for diseased population
  m <- dim(xx)[1]
  # z1 <- xx[, 2]
  # h.z1 <- CVbw(type_kernel = "n", z1, n_pts = 100, seq_bws =NULL)$bw / 2
  pz1 <- rep(0, m)
  for (k in 1:m){
    pz1[k] <- PzIP(xx[k, 2], xx, h.z1) 
  }      
  pz1[which(pz1 < 0.05)] <- 0.05
  # estimate response rate for diseased population
  n <- dim(yy)[1]
  # z2 <- yy[, 2]
  # h.z2 <- CVbw(type_kernel = "n", z2, n_pts = 100, seq_bws =NULL)$bw / 2
  pz2 <- rep(0, n)
  for (k in 1:n){
    pz2[k] <- PzIP(yy[k, 2], yy, h.z2) 
  }
  pz2[which(pz2 < 0.05)] <- 0.05
  a <- ls()
  rm(list=a[which(a != 'pz1' & a != 'pz2')]) 
  list("pz1" = pz1, "pz2" = pz2)
}


# estimated propensity function p(xi) by nonparametrical kernel
PzIP <- function(zz, xx, h.z){
  TINY <- sqrt( .Machine$double.xmin )
  # n <- dim(xx)[1]
  delta <- xx[, 1]
  z <- xx[, 2]
  weight <- dnorm( (z - zz) / h.z)
  if (max(weight) < TINY) weight <- weight + TINY
  weight.all <- sum(weight)
  pz <- sum(weight * delta) / weight.all
  a <- ls()
  rm(list=a[which(a != 'pz')])
  return(pz)  
}

###########                   References                        ################

#Bianco AM, Boente G, González–Manteiga W, Pérez–González A. Estimators for ROC curves with missing biomarkers values and informative covariates. Statistical Methods & Applications. 2023.

#Cheng W, Tang N. Smoothed empirical likelihood inference for ROC curve in the presence of missing biomarker values. Biom J. 2020;62(4):1038-59.

#Long Q, Zhang X, Hsu C-H. Nonparametric multiple imputation for receiver operating characteristics analysis when some biomarker values are missing at random. Stat Med. 2011a;30(26):3149-61.

#Long Q, Zhang X, Johnson BA. Robust estimation of area under ROC curve using auxiliary variables in the presence of missing biomarker values. Biometrics. 2011b;67(2):559-67.

#Wang B, Qin G. Imputation-based empirical likelihood inference for the area under the ROC curve with missing data. Stat Interface. 2012;5(3):319-29.

#Wang B, Qin G. Empirical likelihood-based confidence intervals for the sensitivity of a continuous-scale diagnostic test with missing data. Commun Stat Theory Methods. 2014;43(15):3248-68.
