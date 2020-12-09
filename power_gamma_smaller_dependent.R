library(doParallel)

numCores = 4

registerDoParallel(numCores)

library(EnvStats)
library(actuar)
library(RcppRoll) #for calculating moving sums
library(beepr) #for the signal when the program is finished
library(tidyr)
library(qdapTools)
library(rlist)
library(dplyr)



path_to_sources <- ""


source(paste0(path_to_sources,"MRS.R"))
source(paste0(path_to_sources,"get_sample_dependent.R"))
source(paste0(path_to_sources,"segment_change.R"))
source(paste0(path_to_sources,"MRS_statistic_small.R"))
source(paste0(path_to_sources,"params_list.R"))
source(paste0(path_to_sources,"spots_x_main.R"))

table_crit_val <- read.csv("critical values.csv")

#where to write results
filename='gamma_calc.csv'


#significance level
alpha_ <- c(0.05)

#coefficients of AR or MA process for simulation
AR_par <- c(0.5,0)
MA_par <- c(0,0)


#mean of the structural break
mean_break_list <- c(0.3)

#innovations
innovations <- "pareto"
#hyperparameter of a Pareto distribution: location or Loggamma distr. shapelog
other_par <- 1
#hyperparameter of a Pareto distribution: shape=a or Loggamma distr. ratelog
a <- 2.5
#parameter for symmetrized innovations
p <- 1/2

#size of a 1/4 sample
size_4m <- c(250)

#list of relative lengths of mean change compared to sample size
#l_rel_len <- c(1/30,1/10,1/6,1/4.28, 1/15, 1/7.5,1/5 )
#l_rel_len <- c(1/30,1/10,1/5)
l_rel_len <- c(1/6)

monte_carlo <- 1000

grid_length <- 15
gamma_threshold <- max(0,1/2)
left_side <- seq(from=0, to=gamma_threshold, length.out=grid_length)
right_side <- seq(from=gamma_threshold, to=(gamma_threshold+0.5), by=0.05)
gamma_vector <- left_side[-length(left_side)]



spots_x_main_timed <- function(pars_list){
  system.time(
    results <- spots_x_main(pars_list)
  )
  return(results)
}

parameters <- params_list(AR_par, MA_par, innovations, a, other_par, p, mean_break_list, l_rel_len, gamma)

results <- lapply(parameters, spots_x_main_timed)


long <- list2df(results)
long$X3 <- rep(names(results[[1]]), dim(long)[1]/length(results[[1]]))
wide <- spread(long, X3, X1)

# write csv
if(is.na(file.info(filename)$size)){
  write.csv(wide, 
            filename)
}else{
  write.table(wide, 
              filename, 
              sep = ",", 
              col.names = !file.exists(filename), append = T)
}


