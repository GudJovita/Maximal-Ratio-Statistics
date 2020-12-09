

# paramas
params_list <- function(AR_par, MA_par, innovations, a, other_par, p, mean_break_list, l_rel_len, gamma){
  
  
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

