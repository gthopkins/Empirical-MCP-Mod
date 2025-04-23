# Author: Grant Hopkins
# Date: June 29, 2023
#
# Description: this file makes a function, validate_data(), that confirms the
#              ClinDR conforms to expected format requirements subset at the
#              trial-level. For example, a single trial should have unique
#              doses, identical summary information, and non-missing values for
#              key summary statistics.
#
# Input: a single tibble from md_list in the wrangle_data output
#
# Output: the script provides warnings of unanticipated issues in the data and
#         will return an error to stop the code.


validate_data <- function(data,
                          allow_dup_dose = NULL,
                          allow_differ   = NULL,
                          allow_miss     = NULL,
                          allow_zeros    = NULL) {
  
  if(require(tidyverse)){"You must library the tidyverse package."}
  
  differ <- c("dose", "rslt", "se", "sd", "lcl", "ucl", "alpha", "sampsize",
              "ittsize", "pmiss", "regimen", allow_differ)
  
  not_unique <- data %>%
    select(-matches(differ)) %>%
    map(~ length(unique(.x))) %>%
    bind_rows() %>%
    gather(key = column, value = nunique) %>%
    filter(nunique > 1) %>%
    pull(column)
  
  details <- data %>%
    select(matches(not_unique)) %>%
    map(~ str_c(unique(.x), collapse = ", ")) %>%
    bind_rows() %>%
    gather(key = column, value = levels) %>%
    as.list
  
  miss <- data %>%
    select(dose, rslt, se, sd, sampsize) %>%
    map(~ sum(is.na(.x))) %>%
    bind_rows() %>%
    gather(key = column, value = missing) %>%
    filter(missing > 0) %>%
    pull(column)
  
  zero <- data %>%
    select(rslt, se, sd, sampsize) %>%
    map(~ sum(.x == 0)) %>%
    bind_rows() %>%
    gather(key = column, value = zeros) %>%
    filter(zeros > 0) %>%
    pull(column)
  
  infinite <- data %>%
    summarize(infinite = any(!is.finite(rslt)))
  
  # create an issue flag
  issue <- FALSE
  
  # identify duplicated doses
  if (data %>% pull(dose) %>% duplicated %>% any &
      !(data$trial_id[1] %in% allow_dup_dose)) {
    message("A dose was repeated.")
    issue <- TRUE
  }
  
  # identify columns with unexpectedly non-identical values
  if (length(not_unique) > 0 &
      !(data$trial_id[1] %in% allow_differ)) {
    for (i in 1:length(not_unique)) {
      message("The variable ", details$column[i],
              " unexepectedly has differing values:\n",
              "         >> ", details$levels[i])
      issue <- TRUE
    }
  }
  
  # identify missing values in summary statistics
  if (length(miss) > 0  &
      !(data$trial_id[1] %in% allow_miss)) {
    for (i in 1:length(miss)) {
      message("The variable ", miss[i], " unexepectedly has missing values.")
      issue <- TRUE
    }
  }
  
  # identify zero values in summary statistics
  if (length(zero) > 0 &
      !(data$trial_id[1] %in% allow_zeros)) {
    for (i in 1:length(zero)) {
      message("The variable ", zero[i], " unexepectedly has a zero values.")
      issue <- TRUE
    }
  }
  
  # identify non-finite response
  if (infinite$infinite) {
    message("The response takes non-finite values.")
    issue <- TRUE
  }
  
  if (!issue) {
    cat(paste0("\033[0;", "0;32", "m", "No issues found.", "\033[0m","\n"))
  } else {
    cat(paste0("\033[0;", "0;33", "m", "Address issue.", "\033[0m","\n"))
  }
  
  return(issue)
}
