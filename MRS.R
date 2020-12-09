

#Calculation of maximal ratio statistics

MRS <- function(sample, gamma) {
m <- length(sample)/4

for (u in 1:4){
  #devide sample into 4 equal subsamples
  range_i <- (1+(u-1)*m): (m+(u-1)*m)
  subsample <- sample[range_i]
  
  #assign name to 4 equal subsamples for moving sums as mov.sums.j, where j=1,...,4
  MS.name <- paste("mov.sums.", u, sep = "")
  assign(MS.name, data.frame())  
  
  #create empty data frame for coefficients
  l.coeff <- data.frame()
  
  #assign name to 4 statistics T_nj as T_nj, where j=1,...,4
  T.name <- paste("T_n", u, sep = "")
  assign(T.name, data.frame())  
  
 
  for (l in 1:m) {
    
    #assing name to columns of data frame of mov.sums.l (temporal data frame)  
    nam <- paste("X.centr.", range_i[l], sep = "")
    #calculate moving sums and assign to mov.sums.l 
    assign(nam, RcppRoll::roll_sum(subsample,l,fill=NA, align="left"))
    
    #rename columsn of mov.sums.l
    mov.sums.l <- data.frame(get(nam))
    colnames(mov.sums.l) <- nam
    
    #calculate coefficient
    l.jm <- l^(-gamma)
    
    #merge different rolling windows to same data frame mov.sums.j
    #and coefficients to l.j
    if (l==1){
      
      assign(MS.name, rbind(get(MS.name), mov.sums.l))
      l.coeff<-  rbind(l.coeff, l.jm)
      
    } else {
      
      assign(MS.name, cbind(get(MS.name), mov.sums.l))
      l.coeff<-  cbind(l.coeff, l.jm)
    }
  }
  
  #estimate T_nj, where j=1,...,4
  T_nj <- max(l.coeff)*max(na.omit(as.vector(as.matrix(get(MS.name)))))
  assign(T.name, T_nj)
  
}  
  result <- list(T_n1=T_n1,T_n2=T_n2,T_n3=T_n3,T_n4=T_n4)
  return(result)
  
  
}
