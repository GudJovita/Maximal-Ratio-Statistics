


segment_change <- function(sample, subsample_index_ch, mean_break){
  subsample_change <- sample[subsample_index_ch]+mean_break
  sample[subsample_index_ch] <- subsample_change
  return(sample)
}
