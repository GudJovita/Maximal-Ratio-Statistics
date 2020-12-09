
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
