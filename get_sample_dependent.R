
library(forecast)


get_sample_dependent <- function(AR_par, MA_par, n,  a, other_par, p){
  if(innovations=="pareto"){
    #innovations from Pareto distribution
    X <- EnvStats::rpareto(n, other_par, a )
    X_ <- EnvStats::rpareto(n, other_par, a )
    innov_par <- X-X_
    #sample
    sample <- arima.sim(model=list(ar = c(AR_par[1], AR_par[2]),
                                   ma = c(MA_par[1], MA_par[2])), 
                                   n=n,
                                   innov=innov_par)
  }else if(innovations=="loggamma"){
    X <- actuar::rlgamma(n, shapelog=other_par, ratelog=a)
    #innovations from Uniform distribution
    U <- runif(n, min = 0, max = 1)
    Z <- vector()
    for(i in 1:length(U)){
      if(U[i]<p){
        Z[i] <- 1
      }else{
        Z[i] <- -1
      }
    }
    innov_log <- Z*X
    #sample
    sample <- arima.sim(model=list(ar = c(AR_par[1], AR_par[2]),
                                   ma = c(MA_par[1], MA_par[2])), 
                        n=n,
                        innov=innov_log)
  }
  return(sample)
}
