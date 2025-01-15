## ########################################################################## ##
###################     Simulation Study: simulation        ####################
##                          Author: KS                                        ##
##                      Date: October 2024                                    ##
##                        Out: AUC (ROC)                                      ##
## ########################################################################## ##

###################       Main Simulation                 ######################
rm(list=ls())


###### packages ######
.libPaths("/usw/bbb2836/R-library")
library(mvtnorm)
library(mice)
library(pROC)
library(mitml)
library(mix)
library(doFuture)
library(future)
library(parallelly)



###### paths ######
#path_home <- path_work <- "C:/Users/stahlmann/Documents/Promotion/Simulationsstudie/Analyse/Hummel/Scripte/"
path_work <- "/work/bbb2836/phd/"
path_home <- "/home/bbb2836/"


###### functions ######
source(paste0(path_home,"functions_AUC.R"))


###### set sup simulation ######

# a priori: define missingness pattern (for data generation)
mypattern1 <- matrix(nrow = 1, ncol = 5)
mypattern1[1, ] <- c(0,1,1,1,1) # 1=observed, 0=missing; here missing values in index test (V1)

# number of simulations
nsim <- 1000
nsim 

#define simulation scenarios

# for round 1 ; seed for parallel loop 5273; seed for data generation 4730*i
grid = expand.grid(
  sim = 1:nsim 
  , N = c(500, 1000)
  , p = c(0.1, 0.3, 0.5)
  , AUC_0 = c(0.7, 0.85, 0.9)
  , korr = c(0.5)
  , pm = c(0.1, 0.3, 0.5)
  , mech = c("MCAR", "MAR", "MNAR")
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
writeLines(c("Log of iterations"), paste0(path_work, "log3.txt"))

## simulation: in parallel on multiple nodes
hosts <- system('srun hostname', intern=TRUE)
hosts
hosts_core <- rep(hosts, each=15)
cl <- parallelly::makeClusterPSOCK(hosts_core, rscript_libs = .libPaths())
plan(cluster, workers=cl)

availableCores()
nbrOfWorkers()

suppressMessages({
  
  tfull1 <- proc.time()
  set.seed(5273)
  result <- foreach(i = 1:nrow(grid), .combine = rbind, .options.future=list(seed=TRUE)) %dofuture% {
      
      sink(paste0(path_work, "log3.txt"), append=TRUE)
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
        
        
               
        res <- c(AUC.CCA,CIL.CCA,CIU.CCA,time.CCA,AUC.MI2,CIL.MI2,CIU.MI2,time.MI2,AUC.MIB2,CIL.MIB2,CIU.MIB2,
                 time.MIB2,AUC.HDEL,CIL.HDEL,CIU.HDEL,time.HDEL,AUC.KER,CIL.KER,CIU.KER,time.KER)
  
    }
})

colnames(result) <- c("AUC.CCA","CIL.CCA","CIU.CCA","time.CCA","AUC.MI2","CIL.MI2","CIU.MI2","time.MI2","AUC.MIB2","CIL.MIB2","CIU.MIB2",
                      "time.MIB2","AUC.HDEL","CIL.HDEL","CIU.HDEL","time.HDEL","AUC.KER","CIL.KER","CIU.KER","time.KER")


res <- cbind(grid,result)
tfull2 <- proc.time()
time_sim <- tfull2-tfull1
time_sim

save(res, file = paste0(path_work,"results_part1.Rdata"))

sessionInfo()




