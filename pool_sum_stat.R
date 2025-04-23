# Author: Grant Hopkins
# Date: June 30, 2023
#
# Description: this file pools summary statistics within a trial for different
#              regimens that have a similar total daily dose, as measured by
#              primregimen.
#
# Input: tibble of wrangled metaData grouped by taid, possibly with repeat doses
#
# Output: tibble of wrangled metaData grouped by taid, with any repeat doses
#         pooled into the same summary statistic


pool_sum_stat <- function(df){
  if (nrow(df) == 1) {
    return(df)
  } else{
    # identify the columns to be pooled
    index_sampsize <- which(colnames(df) == "sampsize")
    index_rslt <- which(colnames(df) == "rslt")
    index_sd <- which(colnames(df) == "sd")
    index_se <- which(colnames(df) == "se")
    
    # collect information known from each trial
    all_sampsize <- df[, index_sampsize]
    all_rslt <- df[, index_rslt]
    all_sd <- df[, index_sd]
    
    ################################
    # calculate combined estimates #
    ################################
    
    # combined sample size is the sum of sample size
    pooled_sampsize <- sum(all_sampsize)
    
    # combined sample mean is sum of sample means weighted by sample size
    pooled_rslt <- sum((all_sampsize / sum(all_sampsize)) * all_rslt)
    
    # combined sample variance is calculated using sample variances, sample
    # sizes, and differences of sample means from the combined mean
    pooled_var <- sum((all_sampsize - 1) * all_sd^2 +
                        (all_sampsize) * (all_rslt - pooled_rslt)^2) /
                  (pooled_sampsize - 1)
    
    # combined standard deviation is square root of combined variance
    pooled_sd <- sqrt(pooled_var)
    
    # combined standard error is directly calculated using combined estimates
    pooled_se <- pooled_sd / sqrt(pooled_sampsize)
    
    # insert the pooled estimates
    row <- df[1, ]
    row[, index_sampsize] <- pooled_sampsize
    row[, index_rslt] <- pooled_rslt
    row[, index_sd] <- pooled_sd
    row[, index_se] <- pooled_se
    
    # remove confidence interval data, which was not updated
    row[c("lcl", "ucl", "alpha")] <- NA
    
    return(as.data.frame(row))
  }
}
