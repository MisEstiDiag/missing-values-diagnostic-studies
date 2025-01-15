## ########################################################################## ##
###################     Simulation Study: simulation        ####################
##                          Author: KS                                        ##
##                      Date: October 2023                                    ##
##                        Out: AUC (ROC)                                      ##
## ########################################################################## ##

###### Title: micemix part 2 ######
print("mice lr Part 1")


###################       Main Simulation                 ######################
rm(list=ls())


###### packages ######
#.libPaths("/usw/bbb2836/R-library")
library(mvtnorm)
library(mice)
library(pROC)
library(mitml)
library(doFuture)
library(future)
library(parallelly)



###### paths ######
path_home <- path_work <- "C:/Users/stahlmann/Documents/Promotion/Simulationsstudie/Analyse/hpc/Scripte/"
#path_work <- "/work/bbb2836/test_workdir/"
#path_home <- "/home/bbb2836/"


###### functions ######
source(paste0(path_home,"functions_AUC.R"))


###### set sup simulation ######

# a priori: define missingness pattern (for data generation)
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
writeLines(c("Log of iterations"), paste0(path_work, "log_mice_p.txt"))

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
  plan(multisession, workers=7)
  set.seed(5273)
  result <- foreach(i = 1:nrow(grid), .combine = rbind, .options.future=list(seed=TRUE)) %dofuture% {
      
      sink(paste0(path_work, "log_mice_p.txt"), append=TRUE)
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
      
  
      # MI: mice (>36 Std)
      t1 <- proc.time()
      set.seed(954*i)
      predMatrix <- make.predictorMatrix(data) # in some cases predMatrix should be manually altered
      impMethod <- make.method(data) 
      impMethod[] <- ""
      idx <- names(which(colSums(is.na(data)) > 0))
      impMethod[idx] <- "pmm" 
      invisible(capture.output(imp <- mice(data, m=20, maxit=10, method = impMethod, predictorMatrix = predMatrix))) # imputation
      implist <- mitml::mids2mitml.list(imp)
      fit <- unlist(lapply(implist, function(d) auc(d$D, d$V1))) # estimate AUC
      AUC.mice_p = mean(fit)
      se <- unlist(lapply(implist, function(d) sqrt(var(auc(d$D, d$V1))))) # estimate CI
      var.within = sum(se^2)/length(se) 
      var.between = (sum((fit-AUC.mice_p)^2))/(length(fit)-1)
      se.mice = sqrt(var.within+var.between+(var.between/length(se))) # pooled SE
      CIU.mice_p = AUC.mice_p+1.96*se.mice # confidence interval
      CIL.mice_p = AUC.mice_p-1.96*se.mice
      t2 <- proc.time()
      time_sim <- t2-t1
      time.mice_p <- time_sim[3]
      
      res <- c(AUC.mice_p, CIL.mice_p, CIU.mice_p, time.mice_p)
  
    }
})

colnames(result) <- c("AUC.mice_p","CIL.mice_p","CIU.mice_p","time.mice_p")


res <- cbind(grid,result)
tfull2 <- proc.time()
time_sim <- tfull2-tfull1
time_sim

save(res, file = paste0(path_work,"results_mice_p.Rdata"))

sessionInfo()

