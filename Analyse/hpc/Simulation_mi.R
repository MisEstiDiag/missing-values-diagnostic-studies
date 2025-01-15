## ########################################################################## ##
###################     Simulation Study: simulation        ####################
##                          Author: KS                                        ##
##                      Date: September 2024                                     ##
##                        Out: AUC (ROC)                                      ##
## ########################################################################## ##

###################       Supplement                 ######################
rm(list=ls())


###### packages ######
library(mvtnorm)
library(mice)
library(doFuture)
library(future)
library(parallelly)
library(mi)
library(pROC)


###### functions ######
source("./Analyse/functions_AUC.R")



###### set sup simulation ######

# a priori: define missingness pattern
mypattern1 <- matrix(nrow = 1, ncol = 5)
mypattern1[1, ] <- c(0,1,1,1,1) # 1=observed, 0=missing; here missing values in index test (V1)

# number of simulations
nsim <- 100 # of 100
nsim


#define simulation scenarios
grid = expand.grid(
  sim = 1:nsim 
  , N = c(500)
  , p = c(0.3)
  , AUC_0 = c(0.85)
  , korr = c(0.5)
  , pm = c(0.1, 0.3, 0.5)
  , mech = c("MCAR", "MAR", "MNAR")
)


### AUCmin for power calculation (formula according to Obuchowski et al. 2004
alpha <- qnorm(1-0.05) # alpha set to 5%
beta <- qnorm(1-0.2) # desired power set to 80% / beta=0.1
A <- qnorm(grid$AUC_0)*1.414
V <- (0.0099*exp(-(A^2)/2))*((5*A^2+8)+(A^2+8)/((1-grid$p)/grid$p))
grid$AUC_min <- -(sqrt((alpha*sqrt(0.0792*(1+1/((1-grid$p)/grid$p)))+beta*sqrt(V))^2/(grid$N*grid$p))-grid$AUC_0)

#include a number for the scenario
grid$scenario <- rep(1:(nrow(grid)/nsim), each = nsim)




######  Start Simulation  ######

# document progress: log file
writeLines(c("Log of iterations"),  "./Analyse/log_mi.txt")

suppressMessages({
  
  tfull1 <- proc.time()
  plan(multisession, workers=7)
  set.seed(5273)
  result <-  foreach(i = 1:nrow(grid), .combine = rbind, .options.future=list(seed=TRUE)) %dofuture% {
      
      
      sink("./Analyse/log2.txt", append=TRUE)
      on.exit(sink())
      cat(paste("Starting iteration",i,Sys.time(),Sys.info()[['nodename']],"\n"))
      
        
        # generate data
        # use only data where there is at least 1 observed value for the index test of diseased and non-diseased subjects, respectively
        data <- data.frame(matrix(ncol = 5, nrow = 1))
        colnames(data) <- c("v1", "V2", "V3", "V4", "D")
        data$D <- 0
        seed <- 0
        while (sum(!is.na(data[which(data$D==1), "V1"]))==0 | sum(!is.na(data[which(data$D==0), "V1"]))==0) {
          data = data_fun(N=grid$N[i], p=grid$p[i], AUC_0=grid$AUC_0[i], korr=grid$korr[i], 
                          mech=grid$mech[i], pm=grid$pm[i], mypattern=mypattern1, seed=(4730+seed)*i)
          seed <- seed+1
        }
      
        
        # MI: mi ca 4,8 Tage
        t1 <- proc.time()
        mdf <- missing_data.frame(data)
        mdf <- change(mdf, y = "V1", what = "imputation_method", to = "pmm")
        imp <- mi::mi(mdf, n.chains=20,  n.iter=10, parallel=FALSE, seed=347*i)
        implist <- mi::complete(imp)
        fit <- unlist(lapply(implist, function(d) auc(d$D, d$V1))) # estimate AUC
        AUC.mi_ = mean(fit)
        se <- unlist(lapply(implist, function(d) sqrt(var(auc(d$D, d$V1))))) # estimate CI
        var.within = sum(se^2)/length(se) # 20=no of imputation sets
        var.between = (sum((fit-AUC.mi_)^2))/(length(fit)-1)
        se.mi = sqrt(var.within+var.between+(var.between/length(se))) # pooled SE
        CIU.mi_ = AUC.mi_+1.96*se.mi # confidence interval
        CIL.mi_ = AUC.mi_-1.96*se.mi
        t2 <- proc.time()
        time_sim <- t2-t1
        time.mi_ <- time_sim[3]
        
        
      res <- c(AUC.mi_,CIL.mi_,CIU.mi_,time.mi_)
      
  }

colnames(result) <- c("AUC.mi_","CIL.mi_","CIU.mi_","time.mi_")

})

res <- cbind(grid,result)

tfull2 <- proc.time()
time_sim <- tfull2-tfull1
time_sim

save(res, file = "./Analyse/results_mi.Rdata")

sessionInfo()

