## ########################################################################## ##
###################     Simulation Study: simulation        ####################
##                          Author: KS                                        ##
##                      Date: November 2024                                    ##
##                        Out: AUC (ROC)                                      ##
## ########################################################################## ##

###############         Additional MNAR scenarios             ##################
rm(list=ls())


###### packages ######
#.libPaths("/usw/bbb2836/R-library")
library(mvtnorm)
library(mice)
library(pROC)
library(mitml)
library(mix)
library(doFuture)
library(future)
library(parallelly)



###### paths ######
#path_home <- path_work <- "C:/Users/stahlmann/Documents/Promotion/Simulationsstudie/Analyse/"
#path_work <- "/work/bbb2836/phd/"
#path_home <- "/home/bbb2836/"


###### functions ######
source("./Analyse/functions_AUC.R")

# Description: function to generate simulated data under MNAR, specifying type
# parameters:
# N = sample size
# p = prevalence of target condition
# AUC_0 = true AUC
# korr = correlation between index test and covariates
# mech = missingness mechanism (MCAR, MAR, MNAR)
# pm = proportion of missing values
# mypattern = pattern of missing values (only in the index test)
# seed = seed for data generation

data_mnar <- function(N, p, AUC_0, korr, mech, pm, mypattern, type, seed){
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
  amp <- ampute(data, prop = pm, patterns = mypattern, mech = mech, type = "MID")
  data_mis <- amp$amp
  return(data_mis)
}

###### set sup simulation ######

# a priori: define missingness pattern (for data generation)
mypattern1 <- matrix(nrow = 1, ncol = 5)
mypattern1[1, ] <- c(0,1,1,1,1) # 1=observed, 0=missing; here missing values in index test (V1)

# number of simulations
nsim <- 100 # 1,6 std
nsim 

#define simulation scenarios

# for round 1 ; seed for parallel loop 5273; seed for data generation 4730*i
grid = expand.grid(
  sim = 1:nsim 
  , N = c(500, 1000)
  , p = c(0.1, 0.3, 0.5)
  , AUC_0 = c(0.85)
  , korr = c(0.5)
  , pm = c(0.3)
  , mech = c("MNAR")
  , type = c("MID")
)


### AUCmin for power calculation (formula according to Obuchowski et al. 2004)
alpha <- qnorm(1-0.05) # alpha set to 5%
beta <- qnorm(1-0.2) # desired power set to 80% / beta=0.1
A <- qnorm(grid$AUC_0)*1.414
V <- (0.0099*exp(-(A^2)/2))*((5*A^2+8)+(A^2+8)/((1-grid$p)/grid$p))
grid$AUC_min <- -(sqrt((alpha*sqrt(0.0792*(1+1/((1-grid$p)/grid$p)))+beta*sqrt(V))^2/(grid$N*grid$p))-grid$AUC_0)

#include a number for the scenario
grid$scenario <- rep(1:(nrow(grid)/nsim), each = nsim)




######  Start Simulation  ######

# document progress: log file
writeLines(c("Log of iterations"), paste0(path_work, "log_mnar.txt"))

## simulation: in parallel on multiple nodes
#hosts <- system('srun hostname', intern=TRUE)
#hosts
#hosts_core <- rep(hosts, each=15)
#cl <- parallelly::makeClusterPSOCK(hosts_core, rscript_libs = .libPaths())
#plan(cluster, workers=cl)

availableCores()
nbrOfWorkers()

suppressMessages({
  
  tfull1 <- proc.time()
  plan(multisession, workers=7)
  set.seed(5273)
  result <- foreach(i = 1:nrow(grid), .combine = rbind, .options.future=list(seed=TRUE)) %dofuture% {
      
      sink(paste0(path_work, "log_mnar.txt"), append=TRUE)
      on.exit(sink())
      cat(paste("Starting iteration",i,Sys.time(),Sys.info()[['nodename']],"\n"))
      
      # generate data
      # use only data where there is at least 1 observed value for the index test of diseased and non-diseased subjects, respectively
      data <- data.frame(matrix(ncol = 5, nrow = 1))
      colnames(data) <- c("v1", "V2", "V3", "V4", "D")
      data$D <- 0
      seed <- 0
      while (sum(!is.na(data[which(data$D==1), "V1"]))==0 | sum(!is.na(data[which(data$D==0), "V1"]))==0) {
        data = data_mnar(N=grid$N[i], p=grid$p[i], AUC_0=grid$AUC_0[i], korr=grid$korr[i], 
                        mech=grid$mech[i], pm=grid$pm[i], mypattern=mypattern1, seed=(4730+seed)*i)
        seed <- seed+1
      }
      
  
        # HDEL
        t1 <- proc.time()
        res <- hdel(data=data, d="D", ind="V1", a=0.95, l1=0.5, l2=0.95, tt=seq(0.01, 1-0.01, 0.01), seed=548*i) # warnings indicate that <0.75 there is no logr
        t2 <- proc.time()
        time_sim <- t2-t1
        time.HDEL <- time_sim[3]
        AUC.HDEL <- res$AUC.HDEL1
        CIU.HDEL <- res$CIU.HDEL
        CIL.HDEL <- res$CIL.HDEL
        
        # CCA
        t1 <- proc.time()
        if (sum(!is.na(data[which(data$D==1), "V1"]))==0 | sum(!is.na(data[which(data$D==0), "V1"]))==0) {
          AUC.CCA <- NA # if all subjects with or without target condition have only missing values, this method doesn't work -> set to NA
          CIU.CCA <- NA
          CIL.CCA <- NA
          time.CCA <- NA
        } else {
          ci <- ci.auc(data$D,data$V1)
          t2 <- proc.time()
          time_sim <- t2-t1
          time.CCA <- time_sim[3]
          AUC.CCA <- ci[2]
          CIU.CCA <- ci[3]
          CIL.CCA <- ci[1]
        }
        
        
        # MI2 
        t1 <- proc.time()
        res <- mi2(data=data, d="D", ind="V1", cov_pred=c("V2","V3","V4"), cov_prop=c("V2","V3","V4"), 
                   w=c(0.5,0.5), kN=3, kL=10, tt=seq(0.01, 1-0.01, 0.01), kLevel=0.95, seed=395*i)
        t2 <- proc.time()
        time_sim <- t2-t1
        time.MI2 <- time_sim[3]
        AUC.MI2 <- res$AUC.MI2
        CIU.MI2 <- res$CIU.MI2
        CIL.MI2 <- res$CIL.MI2
        
        
        # KER
        t1 <- proc.time()
        res <- ker(data=data, d="D", ind="V1", cov_prop=c("V2", "V3", "V4"), cov_pred=c("V2", "V3", "V4"), tt=seq(0.01, 1-0.01, 0.01))
        t2 <- proc.time()
        time_sim <- t2-t1
        time.KER <- time_sim[3]
        AUC.KER <- res$AUC.KER
        CIU.KER <- res$CIU.KER
        CIL.KER <- res$CIL.KER
        
        
        # MIB2  
        t1 <- proc.time()
        res <- mi2boot(data=data, d="D", ind="V1", cov_pred=c("V2","V3","V4"), cov_prop=c("V2","V3","V4"), 
                       w=c(0.5,0.5), kN=3, kL=10, tt=seq(0.01, 1-0.01, 0.01), kLevel=0.95, seed=482*i)
        t2 <- proc.time()
        time_sim <- t2-t1
        time.MIB2 <- time_sim[3]
        AUC.MIB2 <- res$AUC.MIB2
        CIU.MIB2 <- res$CIU.MIB2
        CIL.MIB2 <- res$CIL.MIB2
        
        # MI: mice (>36 Std)
        t1 <- proc.time()
        set.seed(954*i)
        predMatrix <- make.predictorMatrix(data) # in some cases predMatrix should be manually altered
        impMethod <- make.method(data) 
        impMethod[] <- ""
        idx <- names(which(colSums(is.na(data)) > 0))
        impMethod[idx] <- "pmm" # pmm for all variables with missing data
        invisible(capture.output(imp <- mice(data, m=20, maxit=10, method = impMethod, predictorMatrix = predMatrix))) # imputation
        implist <- mitml::mids2mitml.list(imp)
        fit <- unlist(lapply(implist, function(d) auc(d$D, d$V1))) # estimate AUC
        AUC.mice = mean(fit)
        se <- unlist(lapply(implist, function(d) sqrt(var(auc(d$D, d$V1))))) # estimate CI
        var.within = sum(se^2)/length(se) 
        var.between = (sum((fit-AUC.mice)^2))/(length(fit)-1)
        se.mice = sqrt(var.within+var.between+(var.between/length(se))) # pooled SE
        CIU.mice = AUC.mice+1.96*se.mice # confidence interval
        CIL.mice = AUC.mice-1.96*se.mice
        t2 <- proc.time()
        time_sim <- t2-t1
        time.mice <- time_sim[3]
        
        
        # MI: mix (>24 Std)
        t1 <- proc.time()
        data_mix <- data[,c(5,2,1,3,4)] # restructure, categorical variables first
        data_mix$D[data_mix$D==0] <- 2 # 2=non-diseased, 1=diseased
        s <- prelim.mix(data_mix, 2) # 1 gives the number of categorical variables in the dataset
        invisible(capture.output(thetahat <- em.mix(s) ))# get starting value
        rngseed(1545*i) # set random number generator seed, must be set (otherwise da.mix wont work), change for simulation
        implist=list()
        for (m in 1:20) { # imputation
          err <- FALSE
          implist[[m]] <- tryCatch({
            invisible(capture.output(theta <- da.mix(s,thetahat,steps=100))) # data augmentation, thetahat provides the starting values
            #print(err)
            implist[[m]] <- as.data.frame(imp.mix(s, theta, data_mix))
          }, 
          error = function(e) {
            err <- TRUE
            #print(err)
            implist[[m]] <- as.data.frame(imp.mix(s, thetahat, data_mix))
          }) 
        }
        fit <- unlist(lapply(implist, function(d) auc(d$D, d$V1))) # estimate AUC
        AUC.mix = mean(fit) 
        se <- unlist(lapply(implist, function(d) sqrt(var(auc(d$D, d$V1))))) # estimate CI
        var.within = sum(se^2)/length(se) # 20=no of imputation sets
        var.between = (sum((fit-AUC.mix)^2))/(length(fit)-1)
        se.mix = sqrt(var.within+var.between+(var.between/length(se))) # pooled SE
        CIU.mix = AUC.mix+1.96*se.mix # confidence interval
        CIL.mix = AUC.mix-1.96*se.mix
        t2 <- proc.time()
        time_sim <- t2-t1
        time.mix <- time_sim[3]
        
        
        # AIPW for defined idx
        t1 <- proc.time()
        res <- aipw(data=data, d="D", ind="V1", cov_prop=c("V2", "V3", "V4"), cov_pred=c("V2", "V3", "V4"), seed=468*i)
        t2 <- proc.time()
        time_sim <- t2-t1
        time.AIPW <- time_sim[3]
        AUC.AIPW <- res$AUC.AIPW
        CIU.AIPW <- res$CIU.AIPW
        CIL.AIPW <- res$CIL.AIPW
        
               
        res <- c(AUC.CCA,CIL.CCA,CIU.CCA,time.CCA,AUC.MI2,CIL.MI2,CIU.MI2,time.MI2,AUC.MIB2,CIL.MIB2,CIU.MIB2,
                 time.MIB2,AUC.HDEL,CIL.HDEL,CIU.HDEL,time.HDEL,AUC.KER,CIL.KER,CIU.KER,time.KER,AUC.mice,CIL.mice,CIU.mice,time.mice,
                 AUC.mix,CIL.mix,CIU.mix,time.mix,AUC.AIPW,CIL.AIPW,CIU.AIPW,time.AIPW)
  
    }
})

colnames(result) <- c("AUC.CCA","CIL.CCA","CIU.CCA","time.CCA","AUC.MI2","CIL.MI2","CIU.MI2","time.MI2","AUC.MIB2","CIL.MIB2","CIU.MIB2",
                      "time.MIB2","AUC.HDEL","CIL.HDEL","CIU.HDEL","time.HDEL","AUC.KER","CIL.KER","CIU.KER","time.KER","AUC.mice","CIL.mice","CIU.mice",
                      "time.mice","AUC.mix","CIL.mix","CIU.mix","time.mix","AUC.AIPW","CIL.AIPW","CIU.AIPW","time.AIPW")


res <- cbind(grid,result)
tfull2 <- proc.time()
time_sim <- tfull2-tfull1
time_sim

save(res, file = "./Analyse/hpc/Simulation_data/results_mnar.Rdata")

sessionInfo()




