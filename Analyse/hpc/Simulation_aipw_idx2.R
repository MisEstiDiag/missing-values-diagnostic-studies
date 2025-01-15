## ########################################################################## ##
###################     Simulation Study: simulation        ####################
##                          Author: KS                                        ##
##                      Date: August 2024                                     ##
##                        Out: AUC (ROC)                                      ##
## ########################################################################## ##

###################       Main Simulation                 ######################
rm(list=ls())


###### packages ######
.libPaths("/usw/bbb2836/R-library")
library(mvtnorm)
library(mice)
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

# a priori: define missingness pattern
mypattern1 <- matrix(nrow = 1, ncol = 5)
mypattern1[1, ] <- c(0,1,1,1,1) # 1=observed, 0=missing; here missing values in index test (V1)

# number of simulations
nsim <- 1000 
nsim 

#define simulation scenarios

# for round2
grid2 = expand.grid(
  sim = 1:nsim 
  , N = c(100, 500, 1000)
  , p = c(0.1, 0.3, 0.5)
  , AUC_0 = c(0.7, 0.85, 0.9)
  , korr = c(0.2, 0.9)
  , pm = c(0.1, 0.3, 0.5)
  , mech = c("MCAR", "MAR", "MNAR")
)

grid3 = expand.grid(
  sim = 1:nsim 
  , N = c(100)
  , p = c(0.1, 0.3, 0.5)
  , AUC_0 = c(0.7, 0.85, 0.9)
  , korr = c(0.5)
  , pm = c(0.1, 0.3, 0.5)
  , mech = c("MCAR", "MAR", "MNAR")
)

grid <- rbind(grid2, grid3) # seed for parallel loop 952; seed for data generation 810*i

### AUCmin for power calculation (formula according to Obuchowski et al. 2004
alpha <- qnorm(1-0.05) # alpha set to 5%
beta <- qnorm(1-0.2) # desired power set to 80% / beta=0.1
A <- qnorm(grid$AUC_0)*1.414
V <- (0.0099*exp(-(A^2)/2))*((5*A^2+8)+(A^2+8)/((1-grid$p)/grid$p))
grid$AUC_min <- -(sqrt((alpha*sqrt(0.0792*(1+1/((1-grid$p)/grid$p)))+beta*sqrt(V))^2/(grid$N*grid$p))-grid$AUC_0)

#include a number for the scenario
grid$scenario <- rep(1:(nrow(grid)/nsim), each = nsim)

# run the aipw function only for defined number of rows (data generation for all rows, so that the seed is correct)
idx.row1 <- round(nrow(grid)/5,0)
idx.row2 <- idx.row1*2
idx.row3 <- idx.row1*3
idx.row4 <- idx.row1*4
idx.row1


######  Start Simulation  ######

# document progress: log file
writeLines(c("Log of iterations"), paste0(path_work, "log.txt"))

## simulation: in parallel on multiple nodes
hosts <- system('srun hostname', intern=TRUE)
hosts
hosts_core <- rep(hosts, each=16)
cl <- parallelly::makeClusterPSOCK(hosts_core, rscript_libs = .libPaths())
plan(cluster, workers=cl)

availableCores()
nbrOfWorkers()

tfull1 <- proc.time()
set.seed(952)
result <-  foreach(i = (idx.row1+1):idx.row2, .combine = rbind, .options.future=list(seed=TRUE)) %dofuture% {
    
    
    sink(paste0(path_work, "log.txt"), append=TRUE)
    on.exit(sink())
    cat(paste("Starting iteration",i,Sys.info()[['nodename']],"\n"))
    
      
      # generate data
      # use only data where there is at least 1 observed value for the index test of diseased and non-diseased subjects, respectively
      data <- data.frame(matrix(ncol = 5, nrow = 1))
      colnames(data) <- c("v1", "V2", "V3", "V4", "D")
      data$D <- 0
      seed <- 0
      while (sum(!is.na(data[which(data$D==1), "V1"]))==0 | sum(!is.na(data[which(data$D==0), "V1"]))==0) {
        data = data_fun(N=grid$N[i], p=grid$p[i], AUC_0=grid$AUC_0[i], korr=grid$korr[i], 
                        mech=grid$mech[i], pm=grid$pm[i], mypattern=mypattern1, seed=(810+seed)*i)
        seed <- seed+1
      }
    
      
      # AIPW for defined idx
      t1 <- proc.time()
      res <- aipw(data=data, d="D", ind="V1", cov_prop=c("V2", "V3", "V4"), cov_pred=c("V2", "V3", "V4"), seed=468*i)
      t2 <- proc.time()
      time_sim <- t2-t1
      time.AIPW <- time_sim[3]
      AUC.AIPW <- res$AUC.AIPW
      CIU.AIPW <- res$CIU.AIPW
      CIL.AIPW <- res$CIL.AIPW
      
    res <- c(i,AUC.AIPW,CIL.AIPW,CIU.AIPW,time.AIPW)
    
}

colnames(result) <- c("Row.grid","AUC.AIPW","CIL.AIPW","CIU.AIPW","time.AIPW")

res <- cbind(grid[(idx.row1+1):idx.row2,],result)

tfull2 <- proc.time()
time_sim <- tfull2-tfull1
time_sim

save(res, file = paste0(path_work,"results_aipw2.Rdata"))

sessionInfo()

