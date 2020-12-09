

MRS_statistic_small <- function(sample, gamma, alpha, X_crit_1){

          
          MRS_res <- MRS(sample, gamma)
          RMS <- max(MRS_res$T_n1/MRS_res$T_n3, MRS_res$T_n3/MRS_res$T_n1, MRS_res$T_n2/MRS_res$T_n4, MRS_res$T_n4/MRS_res$T_n2)
          
          test_result <- ifelse(RMS>X_crit_1,1,0)
          return(test_result)


}