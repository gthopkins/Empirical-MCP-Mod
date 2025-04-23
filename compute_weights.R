# Author: Charles Liu; small adjustments by Grant Hopkins
# Date: July 10, 2023
#
# Description: the script produces weights of the information criteria to be
#              used in subsequent analysis and prediction.
#
# Input: vector of information criteria from different model assumptions
#
# Output: model weights in same order as provided

compute_weights <- function(IC, keep_name = TRUE){ 
  delta_IC <- IC - min(IC)
  weights <- exp(-0.5 * delta_IC) / sum(exp(-0.5 * delta_IC))
  if (keep_name) {
    names(weights) <- names(IC)
  }
  return(weights)
}