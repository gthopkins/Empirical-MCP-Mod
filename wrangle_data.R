# Author: Grant Hopkins
# Date: August 14, 2023
#
# Description: this file wrangles the data included in the ClinDR package such
#              that the studies are grouped by drug within a therapeutic area
#              within a protocol (e.g., a specific clinical trial). Specific
#              data subsetting and manipulations are described in comments.
#
# Input: metaData file from clinDR package
#
# Output: the script returns a tibble containing the wrangled data


wrangle_data <- function(data, min_num_dose = 5){
  if(require(tidyverse)){"You must library the tidyverse package."}
  
  # read in the data
  data <- tibble(data)
  
  # process the data
  md_proc <- data %>%
    
    # identify each unique drug, therapeutic area, and protocol combination
    mutate(label = "id") %>%
    unite(temp_id, c("taid", "protid"), sep = "_", remove = FALSE) %>%
    unite(trial_id, c("label", "temp_id"), sep = "", remove = TRUE) %>%
    
    # keep only small molecule drugs
    filter(drugtype == "SMALL MOLECULE") %>%
    
    # consider FDA 2009-2014 & 2014-2017, and Pfizer 1998-2009 & 2009-2018
    filter(metasource %in% c("FDA1417", "FDA914", "PF09", "PFIZERUPDATE18")) %>%
    
    # keep only the primary endpoint
    filter(etype == 1) %>%
    
    # keep only the primary regimen
    filter(primregimen == 1) %>%
    
    # keep only binary and continuous endpoints (not time-to-event)
    filter(primtype %in% c("BINARY", "CONTINUOUS")) %>%
    
    # keep only sub-population evaluated on the most doses
    filter(poptype == 1) %>%
    
    # pool response means together when total daily dose (TDD) is equivalent
    group_by(trial_id, dose) %>%
      group_modify(~ {.x %>% pool_sum_stat}) %>%
    ungroup() %>%
    
    # now assess each trial individually
    group_by(trial_id) %>%  
    
      # keep only trials that have an active placebo group
      filter(any(dose == 0)) %>%
      
      # identify number of unique dose-levels for each trials
      mutate(ndose = length(unique(dose))) %>%
    
    # stop grouping by trial_id
    ungroup() %>%
    
    # keep trials with specified minimum number of doses (5 by default)
    filter(ndose >= min_num_dose)
  
  # split data into list of tibbles by trial_id
  md_list <- split(md_proc, md_proc$trial_id)
  
  return(list(md_list = md_list, md_proc = md_proc))
}
