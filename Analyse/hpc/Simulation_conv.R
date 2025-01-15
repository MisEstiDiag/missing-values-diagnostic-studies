## ########################################################################## ##
###################     Simulation Study: simulation        ####################
##                          Author: KS                                        ##
##                      Date: 28 Apr 2023                                     ##
##                        Out: AUC (ROC)                                      ##
## ########################################################################## ##


###### Title: CONV ######
print("CONV")


###################       Main Simulation                 ######################
rm(list=ls())

###### packages ######
.libPaths("/usw/bbb2836/R-library")
library(mvtnorm)
library(mice)
library(doFuture)
library(future)
library(parallelly)
library(pROC)



###### paths ######
#path_home <- path_work <- "C:/Users/stahlmann/Documents/Promotion/Simulationsstudie/Analyse/Hummel/Scripte/"
path_work <- "/work/bbb2836/phd/"
path_home <- "/home/bbb2836/"


###### functions ######
source(paste0(path_home,"functions_AUC.R"))


###### set sup simulation ######
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

### AUCmin for power calculation (formula according to Obuchowski et al. 2004)
alpha <- qnorm(1-0.05) # alpha set to 5%
beta <- qnorm(1-0.2) # desired power set to 80% / beta=0.2
A <- qnorm(grid$AUC_0)*1.414
V <- (0.0099*exp(-(A^2)/2))*((5*A^2+8)+(A^2+8)/((1-grid$p)/grid$p))
grid$AUC_min <- -(sqrt((alpha*sqrt(0.0792*(1+1/((1-grid$p)/grid$p)))+beta*sqrt(V))^2/(grid$N*grid$p))-grid$AUC_0)

#include a number for the scenario
grid$scenario <- rep(1:(nrow(grid)/nsim), each = nsim)


######  Start Simulation  ######

# document progress: log file
writeLines(c("Log of iterations"), paste0(path_work, "log_conv.txt"))

  
  ## simulation: in parallel on multiple nodes
  hosts <- system('srun hostname', intern=TRUE)
  hosts
  hosts_core <- rep(hosts, each=16)
  cl <- parallelly::makeClusterPSOCK(hosts_core, rscript_libs = .libPaths())
  plan(cluster, workers=cl)
  
  availableCores()
  nbrOfWorkers()

suppressMessages({  

  tfull1 <- proc.time()

  set.seed(5273)
  result <- foreach(i = 1:nrow(grid), .combine = rbind, .options.future=list(seed=TRUE)) %dofuture% {
    
    sink(paste0(path_work, "log_ipl.txt"), append=TRUE)
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
      
      # CONV
      t1 <- proc.time()
      res <- conv(data=data, d="D", ind="V1", cov_prop=c("V2", "V3", "V4"), cov_pred=c("V2", "V3", "V4"), tt=seq(0.01, 1-0.01, 0.01))
      t2 <- proc.time()
      time_sim <- t2-t1
      time.CONV <- time_sim[3]
      AUC.CONV <- res$AUC.CONV
      CIU.CONV <- res$CIU.CONV
      CIL.CONV <- res$CIL.CONV
      
  
    res <- c(AUC.CCA,CIL.CCA,CIU.CCA,time.CCA, AUC.CONV,CIL.CONV,CIU.CONV,time.CONV)
             
  }

  colnames(result) <- c("AUC.CCA","CIL.CCA","CIU.CCA","time.CCA","AUC.CONV","CIL.CONV","CIU.CONV","time.CONV")
                        
})

res <- cbind(grid,result)
tfull2 <- proc.time()
time_sim <- tfull2-tfull1
time_sim

save(res, file = paste0(path_work,"results_conv.Rdata"))

sessionInfo()





