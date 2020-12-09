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

path_to_sources <- "C:/Users/JG/Documents/PhD1/Straipsniai/Mean changes by MRS/Project/dependent case/"

source(paste0(path_to_sources,"MRS.R"))
source(paste0(path_to_sources,"get_sample_dependent.R"))
source(paste0(path_to_sources,"segment_change.R"))
source(paste0(path_to_sources,"MRS_statistic_small.R"))


table_crit_val <- read.csv("C:/Users/JG/Documents/PhD1/Straipsniai/Mean changes by MRS/Project/critical values.csv")

#where to write results
#filename='C:/Users/JG/Documents/PhD/Straipsniai/Mean changes by MRS/Results/v2/results_power_paral_small.csv'
#filename='C:/Users/JG/Documents/PhD/Straipsniai/Mean changes by MRS/Results/v2/results_power_table_test.csv'
filename='C:/Users/JG/Documents/PhD1/Straipsniai/Mean changes by MRS/Results/dependent_case/Figure1/calc_1.csv'



#significance level
alpha_ <- c(0.05)

#hyperparameter of a Pareto distribution: location or Loggamma distr. shapelog
other_par <- 1
#hyperparameter of a Pareto distribution: shape=a or Loggamma distr. ratelog
a <- 2.5


#coefficients of AR or MA process for simulation
AR_par <- c(0.5,0)
MA_par <- c(0,0)
sd_par <- 1
#mean of the structural break
mean_break_list <- c(0.2, 0.4, 0.6, 0.8,1,1.2,1.4,1.6,1.8,2)

#innovations
innovations <- "pareto"

#parameter for symmetrized innovations
p <- 1/2
#gamma > max(0,1/2-1/a)
#gamma_threshold <- max(0,1/2-1/a)
#gamma <- round(gamma_threshold-0.01,2)
#gamma <- 0.1
#size of a 1/4 sample
#size_4m <- c(75)
size_4m <- c(75)
#size_4m <- c(150)
#list of relative lengths of mean change compared to sample size
l_rel_len <- c(1/30,1/15, 1/10, 1/7.5,1/6, 1/5 ,1/4.28)
monte_carlo <- 1000

#grid_length <- 10
#gamma_threshold <- max(0,1/2-1/a)
#left_side <- seq(from=0, to=gamma_threshold, length.out=grid_length)
#right_side <- seq(from=gamma_threshold, to=(gamma_threshold+0.5), by=0.05)
#gamma_vector <- left_side[-length(left_side)]

# gamma_vector <- c(round(gamma_threshold-0.01,2), round(gamma_threshold-0.1,2), 
#                   round(gamma_threshold-0.2,2), round(gamma_threshold-0.29,2))
gamma_vector <- c(0.09)



# paramas
params_list <- function(AR_par, MA_par, innovations, a, other_par, p, mean_break_list, l_rel_len, gamma){
  
  # lappend <- function (lst, ...){
  #   lst <- c(lst, list(...))
  #   return(lst)
  # }
  # defined_params <- list()
  # defined_params <- c(distribution, other_par, a, p, mean_break, l_rel_len, gamma)
  params <- list() # alpha, m, n, k
  
  i <- 1
  
  for(g in 1:length(gamma_vector)){ 
    gamma <- gamma_vector[g]
  
  for(c in 1:length(alpha_)){
    alpha <- alpha_[c]
    
    
    for(mean_break in mean_break_list){
      #j denotes a sample of different size
      for (j in 1:length(size_4m)){
        
        m=size_4m[j]
        #sample size
        n <- 4*m
        
        for (k in 1:length(l_rel_len)){
          # temp_params <- c(alpha, m, n, k)
          # temp_params <- c(temp_params, defined_params)
          # params <- lappend(params, temp_params)
          
          temp_params <- list()
          temp_params$alpha <- alpha
          temp_params$m <- m
          temp_params$n <- n
          temp_params$k <- k
          temp_params$AR_par <- AR_par
          temp_params$MA_par <- MA_par
         
          temp_params$innovations <- innovations
          temp_params$a <- a
          temp_params$other_par <- other_par
          temp_params$p <- p
          
          temp_params$mean_break <- mean_break
          temp_params$l_rel_len <- l_rel_len
          temp_params$gamma <- gamma
          
          
          params[[i]] <- temp_params
          i <- i+1
        }
      }
    }
  }
  
  }
  return(params)
}




# calculation
spots_x_main <- function(pars_list){
  
  alpha <- pars_list$alpha
  m <- pars_list$m
  n <- pars_list$n
  k <- pars_list$k
  

  AR_par <- pars_list$AR_par
  MA_par <- pars_list$MA_par
  
  innovations <- pars_list$innovations
  a <- pars_list$a
  other_par <- pars_list$other_par
  p <- pars_list$p
  
  mean_break <- pars_list$mean_break
  l_rel_len <- pars_list$l_rel_len
  gamma <- pars_list$gamma
  
  
  
  message(pars_list)

  alpha_filter <- alpha
  gamma_filter <- round(gamma*2,1)/2
  
  X_crit_1 <- table_crit_val %>% 
    filter(alpha %in% alpha_filter & gamma %in% gamma_filter ) %>%
    select(critical_val) 
  
  event_space_list <- list()
  window <- floor(n*l_rel_len[k])-1
  
  monte_i <- 0
  while (monte_i < monte_carlo) {
    message(monte_i)
  for(s in 1:4){
    
    spots <- seq(from=((s-1)*m+1), to=(m*s-window), by=1)
    
    
    # moving window
    for(i in 1:length(spots)){
      
      
      # main function
      # test_result
      
      subsample_index_ch <- spots[i]:(spots[i]+window)
      
      #subsample_ind <- seq(1:m)[-subsample_index_ch]      
      
      ##########################################################
      sample <- get_sample_dependent(AR_par, MA_par, n, a, other_par, p)
      ############################################################
      sample <- segment_change(sample, subsample_index_ch, mean_break)
      
      # plot_sample(sample, n, l_rel_len, subsample_index_ch, alpha, k)
    
      test_result <- MRS_statistic_small(sample, gamma, alpha, X_crit_1)
      
      #ifelse(RMS>X_crit_2,1,0)
      event_space_list <- append(event_space_list, test_result)
      monte_i <- monte_i+1
      
    }	
  }
  
  }
  
  
  results <- list()
  
  results$event_space <- Reduce("+",event_space_list)
  results$sample_space <- list.count(event_space_list)
  results$power <- results$event_space/results$sample_space
  results$l <- floor(n*l_rel_len[k])
  results$alpha <- alpha
  results$m <- m
  results$n <- n
  
  results$AR1 <- AR_par[1]
  results$AR2 <- AR_par[2]
  results$MA1 <- MA_par[1]
  results$MA2 <- MA_par[2]
  results$innovations <- innovations
  results$a <- a
  results$other_par <- other_par
  results$p <- p
  
  results$mean_break <- mean_break
  results$gamma <- gamma

  
  return(results)
}

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

