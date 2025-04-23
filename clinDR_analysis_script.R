# Author: Grant Hopkins
# Date: January 19, 2024
#
# Description: this script produces the figures included in the manuscript, with
#              additional sensitivity analysis 
#
# Input: all functions included in the GitHub folder
#
# Output: plots produced in folders that will be constructed and additional checks
#         made in the R console

####################
# Library packages #
####################

# general analysis package
library(tidyverse) # data wrangling

# MCP-Mod packages
library(clinDR) # gather metaData
library(DoseFinding) # implement MCP-Mod curves

# plot packages
library(pheatmap) # hierarchical clustered heatmap
library(ggh4x) # scale_y_dendrogram on ggplot
library(dendextend) # flip y-axis of dendrogram
library(ggplotify) # convert pheatmap to ggplot

# table packages
library(matrixStats) # produce row ranks


####################
# Attach function #
####################

# general functions
source("pool_sum_stat.R")
source("wrangle_data.R")
source("validate_data.R")
source("compute_weights.R")

# functions for continuous endpoints
source("make_fits_cont.R")
source("make_plot_cont.R")

# functions for binary endpoints
source("make_fits_bin.R")
source("make_plot_bin.R")

#############################################
# Fit models under two settings:            #
#   1) a minimum of four doses in the study #
#   2) a minimum of five doses in the study #
#############################################

for (ndose in 4:5) {
  dosage <- case_when(ndose == 4 ~ "4dose",
                      ndose == 5 ~ "5dose")
  
  ################
  # Wrangle data #
  ################
  
  data(metaData)
  md_list <- wrangle_data(metaData, min_num_dose = ndose)$md_list
  md_proc <- wrangle_data(metaData, min_num_dose = ndose)$md_proc
  rm(list = c("metaData"))
  
  ################################
  # Manually correct/adjust data #
  ################################
  
  # manually edit id1012_1 with non-placebo-adjusted raw data found via
  # https://www.tga.gov.au/sites/default/files/auspar-macitentan-140428-cer.docx
  
  md_list[["id1012_1"]] <- md_list[["id1012_1"]] %>%
    arrange(dose) %>%
    mutate(dose     = c(  0.0,   0.3,    1.0,    3.0,   10.0)) %>%
    mutate(rslt     = c( -7.9,  -8.8,  -10.0,  -10.8,  -11.8)) %>%
    mutate(sd       = c(  8.0,   8.6,    6.7,    7.3,    6.5)) %>%
    mutate(se       = c(  1.1,   1.2,    0.9,    1.0,    0.9)) %>%
    mutate(sampsize = c( 54  ,  54  ,   60  ,   57  ,   56  )) %>%
    mutate(lcl      = c(-10.1, -11.2,  -11.7,  -12.7,  -13.5)) %>%
    mutate(ucl      = c( -5.7,  -6.4,   -8.3,   -8.8,  -10.0))
  
  md_proc[md_proc$trial_id == "id1012_1", ] <- md_list[["id1012_1"]]
  
  # manually recalculate the true sample proportion, lost to rounding
  
  if (all(c("id1001_1", "id1064_2") %in% names(md_list))) {
    md_list[["id1001_1"]] <- md_list[["id1001_1"]] %>%
      rowwise() %>%
      mutate(rslt = ifelse(rslt <= 0.5,
                           min(Re(polyroot(c(-se^2 * sampsize, 1, -1)))),
                           max(Re(polyroot(c(-se^2 * sampsize, 1, -1)))))) %>%
      ungroup()
    
    md_proc[md_proc$trial_id == "id1001_1", ] <- md_list[["id1001_1"]]
    
    md_list[["id1064_2"]] <- md_list[["id1064_2"]] %>%
      rowwise() %>%
      mutate(rslt = ifelse(rslt <= 0.5,
                           min(Re(polyroot(c(-se^2 * sampsize, 1, -1)))),
                           max(Re(polyroot(c(-se^2 * sampsize, 1, -1)))))) %>%
      ungroup()
    
    md_proc[md_proc$trial_id == "id1064_2", ] <- md_list[["id1064_2"]]
  }
  
  ##############################################
  # Fit models under two settings:             #
  #   1) default exponential parameter space   #
  #   2) extended exponential  parameter space #
  ##############################################
  
  for (exp_set in c("defexp", "extexp")) {
    def_exp <- case_when(exp_set == "defexp" ~ TRUE,
                         exp_set == "extexp" ~ FALSE)
    
    # make data frame to store results
    results <- matrix(NA, nrow = length(md_list), ncol = 24)
    colnames(results) <- c("taid", "protid", "sum_sampsize",
                           "eMax_maxLL",      "eMax_AIC",      "eMax_BIC",
                           "linLog_maxLL",    "linLog_AIC",    "linLog_BIC",
                           "quad_maxLL",      "quad_AIC",      "quad_BIC",
                           "lin_maxLL",       "lin_AIC",       "lin_BIC",
                           "exp_maxLL",       "exp_AIC",       "exp_BIC",
                           "intercept_maxLL", "intercept_AIC", "intercept_BIC",
                           "sigEmax_maxLL",   "sigEmax_AIC",   "sigEmax_BIC")
    
    # iterate through the list of datasets
    for (i in 1:length(md_list)) {
      
      # print the trial id (in bold)
      cat(paste0("\033[0;", "1;39", "m",
                 i, ". Dataset ", names(md_list)[i],
                 "\033[0m","\n"))
      
      # subset to the dataset
      data <- md_list[[i]]
      
      # perform data validation
      flag <- validate_data(data, allow_zeros = c("id23_3", "id1001_1", "id1064_2"))
      if (flag) {next}
      
      # prepare sub-directories to save main plots and sensitivity analyses
      dir.create(file.path(getwd(), file.path("plots")), showWarnings = FALSE)
      dir.create(file.path(getwd(), file.path("plots", exp_set)), showWarnings = FALSE)
      dir.create(file.path(getwd(), file.path("plots", exp_set, dosage)), showWarnings = FALSE)
      
      # fit model
      if (data$primtype[1] == "BINARY") {
        fits <- make_fits_bin(data, default_exp = def_exp)
        
        dir.create(file.path(getwd(), file.path("plots", exp_set, dosage, "bin")),
                   showWarnings = FALSE)
        tiff(filename = paste0(file.path("./plots", exp_set, dosage, "bin"),
                               "/plot_id", data$taid[1], "_", data$protid[1], ".tiff"),
             height = 8, width = 11, units = 'in', res = 300, compression = "lzw")
        print(make_plot_bin(fits = fits, data = data)$plot)
        dev.off()
        
      } else if (data$primtype[1] == "CONTINUOUS") {
        fits <- make_fits_cont(data, default_exp = def_exp)
        
        dir.create(file.path(getwd(), file.path("plots", exp_set, dosage, "cont")),
                   showWarnings = FALSE)
        tiff(filename = paste0(file.path("./plots", exp_set, dosage, "cont"),
                               "/plot_id", data$taid[1], "_", data$protid[1], ".tiff"),
             height = 8, width = 11, units = 'in', res = 300, compression = "lzw")
        print(make_plot_cont(fits = fits, data = data)$plot)
        dev.off()
        
      } else {
        message("The endpoint is neither binary nor continuous.")
      }
      
      
      # save results
      results[i, "taid"] <- data$taid[!is.na(data$taid)][1]
      results[i, "protid"] <- data$protid[!is.na(data$protid)][1]
      results[i, "sum_sampsize"] <- sum(data$sampsize)
      
      results[i, "eMax_maxLL"] <- fits$eMax$maxLL
      results[i, "eMax_AIC"] <- fits$eMax$AIC
      results[i, "eMax_BIC"] <- fits$eMax$BIC
      
      results[i, "linLog_maxLL"] <- fits$linLog$maxLL
      results[i, "linLog_AIC"] <- fits$linLog$AIC
      results[i, "linLog_BIC"] <- fits$linLog$BIC
      
      results[i, "quad_maxLL"] <- fits$quad$maxLL
      results[i, "quad_AIC"] <- fits$quad$AIC
      results[i, "quad_BIC"] <- fits$quad$BIC
      
      results[i, "lin_maxLL"] <- fits$lin$maxLL
      results[i, "lin_AIC"] <- fits$lin$AIC
      results[i, "lin_BIC"] <- fits$lin$BIC
      
      results[i, "exp_maxLL"] <- fits$exp$maxLL
      results[i, "exp_AIC"] <- fits$exp$AIC
      results[i, "exp_BIC"] <- fits$exp$BIC
      
      results[i, "intercept_maxLL"] <- fits$intercept$maxLL
      results[i, "intercept_AIC"] <- fits$intercept$AIC
      results[i, "intercept_BIC"] <- fits$intercept$BIC
      
      results[i, "sigEmax_maxLL"] <- fits$sigEmax$maxLL
      results[i, "sigEmax_AIC"] <- fits$sigEmax$AIC
      results[i, "sigEmax_BIC"] <- fits$sigEmax$BIC
      
      # remove temporary objects
      rm(list = c("data", "fits", "flag"))
    }
    
    ############################
    # Clean the results output #
    ############################
    
    # summarize the results
    results_adj <- as_tibble(results) %>%
      
      # reconstruct the ID
      mutate(id = "id") %>%
      unite(trial_id, c(id, taid), sep = "", remove = FALSE) %>%
      unite(trial_id, c(trial_id, protid), remove = FALSE) %>%
      
      # calculate difference in AIC between Emax and other models
      mutate(del_AIC_linLog  = eMax_AIC - linLog_AIC,
             del_AIC_lin     = eMax_AIC - lin_AIC,
             del_AIC_quad    = eMax_AIC - quad_AIC,
             del_AIC_exp     = eMax_AIC - exp_AIC,
             del_AIC_sigEmax = eMax_AIC - sigEmax_AIC) %>%
      
      # identify when Emax performs better or worse pairwise
      mutate(linLog_better  = ifelse(del_AIC_linLog > 0, "Log-linear", "Emax"),
             quad_better    = ifelse(del_AIC_quad > 0, "Quadratic", "Emax"),
             lin_better     = ifelse(del_AIC_lin > 0, "Linear", "Emax"),
             exp_better     = ifelse(del_AIC_exp > 0, "Exponential", "Emax"),
             sigEmax_better = ifelse(del_AIC_sigEmax > 0, "Sigmoid-Emax",
                                     "Hyperbolic-Emax")) %>%
      
      # in each row...
      rowwise() %>%
      
      # determine which model has the lowest AIC
      mutate(best_AIC_mod = names(which.min(pick(eMax_AIC, linLog_AIC, quad_AIC,
                                                 lin_AIC, exp_AIC)))) %>%
      
      # retrieve the lowest AIC among all models, i.e. "best AIC"
      mutate(best_AIC = min(pick(eMax_AIC, linLog_AIC, quad_AIC,
                                 lin_AIC, exp_AIC))) %>%
      
      # calculate the difference in the null AIC and the "best AIC"
      mutate(del_minAIC_intercept = best_AIC - intercept_AIC) %>%
      
      # make a threshold such that we filter models where null AIC is comparable
      mutate(null_win = ifelse(del_minAIC_intercept >= -2.2, TRUE, FALSE)) %>%
      
      # now stop grouping by row
      ungroup()
    
    ##########################
    # Remove strong null fit #
    ##########################
    
    # identify which trials have comparable intercept fit
    results_adj %>%
      filter(null_win) %>%
      arrange(desc(del_minAIC_intercept)) %>%
      select(trial_id, best_AIC, intercept_AIC, del_minAIC_intercept)
    
    # remove trials with comparable intercept fit
    results_adj <- results_adj %>%
      filter(!null_win)
    
    ##########################
    # Remove duplicated data #
    ##########################
    
    # data is identical to id1007_1
    results_adj <- results_adj %>%
      filter(trial_id != "id1015_1")
    
    #########################################
    # Produce plots for default exponential #
    #########################################
    
    # establish plot colors
    my_cols <- c("#000000", "#56B4E9", "#CC79A7", "#F0E442",
                 "#009E73", "#D55E00", "#0072B2", "#E69F00")
    
    #########################
    # Calculate AIC weights #
    #########################
    
    # AIC weights
    weights_AIC <- data.frame(matrix(nrow = nrow(results_adj), ncol = 7))
    colnames(weights_AIC) <- c("trial_id", "sum_sampsize",
                               "eMax_weight", "linLog_weight",
                               "quad_weight", "lin_weight", 
                               "exp_weight")
    for (i in 1:nrow(results_adj)) {
      weights_AIC[i, 1:2] <- results_adj[i, c("trial_id", "sum_sampsize")]
      weights_AIC[i, 3:7] <- compute_weights(results_adj[i, c("eMax_AIC", "linLog_AIC",
                                                              "quad_AIC", "lin_AIC", 
                                                              "exp_AIC")])
    }
    
    # make numerical weight matrix
    weights_matrix <- weights_AIC[, 3:7]
    rownames(weights_matrix) <- weights_AIC$trial_id
    colnames(weights_matrix) <- c("Emax", "Log-Linear",
                                  "Quadratic", "Linear", 
                                  "Exponential")
    
    # wrangle weights into long format
    weights_AIC <- weights_AIC %>%
      arrange(desc(eMax_weight)) %>%
      mutate(trial_id = fct_reorder(trial_id, eMax_weight)) %>%
      pivot_longer(contains("_weight"), names_to = "Model", values_to = "AIC_weight") %>%
      mutate(Model = as.factor(as.character(Model))) %>%
      mutate(Model = recode_factor(Model,
                                   'eMax_weight'   = 'Emax',
                                   'linLog_weight' = 'Log-linear',
                                   'quad_weight'   = 'Quadratic',
                                   'lin_weight'    = 'Linear',
                                   'exp_weight'    = 'Exponential'))
    
    ##################
    # AIC Rank Table #
    ##################
    
    tab <- results_adj %>%
      as_tibble %>%
      select(eMax_AIC, linLog_AIC, lin_AIC, quad_AIC, exp_AIC) %>%
      as.matrix %>%
      rowRanks(x = ., ties.method = "max") %>%
      as_tibble %>%
      pivot_longer(everything(), names_to = "Model", values_to = "Rank") %>%
      mutate(Model = recode_factor(Model,
                                   'eMax_AIC' = 'Emax',
                                   'linLog_AIC' = 'Log-linear',
                                   'lin_AIC' = 'Linear',
                                   'quad_AIC' = 'Quadratic',
                                   'exp_AIC' = "Exponential")) %>%
      dplyr::count(Model, Rank) %>%
      pivot_wider(names_from = Rank, names_prefix = "Rank ",
                  values_from = n, values_fill = 0) %>%
      mutate(
        across(starts_with("Rank"), ~ ifelse(.x == sort(.x, partial = 5)[5], 0, 1), .names = "freq5_{col}"),
        across(starts_with("Rank"), ~ ifelse(.x == sort(.x, partial = 4)[4], 0, 1), .names = "freq4_{col}"),
        across(starts_with("Rank"), ~ ifelse(.x == sort(.x, partial = 3)[3], 0, 1), .names = "freq3_{col}"),
        across(starts_with("Rank"), ~ ifelse(.x == sort(.x, partial = 2)[2], 0, 1), .names = "freq2_{col}"),
        across(starts_with("Rank"), ~ ifelse(.x == sort(.x, partial = 1)[1], 0, 1), .names = "freq1_{col}")
      ) %>%
      arrange(pick(starts_with("freq"))) %>%
      select(Model, starts_with("Rank")) %>%
      mutate(y = 1:n()) %>%
      pivot_longer(cols = starts_with("Rank"), names_to = "Rank", values_to = "value") %>%
      mutate(x = -as.numeric(as.character(gsub("Rank ", "", Rank)))) %>%
      mutate(Rank = fct_reorder(Rank, -x),
             Model = fct_reorder(Model, -y))
    
    tiff(filename = file.path("./plots", exp_set, dosage, "model_ranks_AIC.tiff"),
         height = 2, width = 4, units = 'in', res = 300, compression = "lzw")
    print(ggplot(tab, aes(x = Rank, y = Model)) +
            geom_tile(aes(fill = value)) +
            geom_text(aes(label = value)) +
            scale_fill_gradient2(low = "white", 
                                 mid = "white", 
                                 high = "#228B22", 
                                 midpoint = 0) + # determine the colour
            theme(panel.grid.major.x = element_blank(), #no gridlines
                  panel.grid.minor.x = element_blank(), 
                  panel.grid.major.y = element_blank(), 
                  panel.grid.minor.y = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.ticks.y = element_blank(),
                  panel.background=element_rect(fill = "white"),
                  axis.text.x = element_text(size = 8, face = "bold"),
                  axis.text.y = element_text(size = 8, face = "bold")) + 
            theme(legend.title = element_text(face = "bold", size = 14),
                  legend.position = "none") + 
            scale_x_discrete(name = "", position = "top") +
            scale_y_discrete(name = ""))
    dev.off()
    
    ###############
    # AIC Heatmap #
    ###############
    
    # cluster rows under Hellinger distance where weights are considered PMFs
    clust_rows <- hclust(dist(sqrt(weights_matrix), method = "euclidean") / sqrt(2),
                         method = "complete")
    
    # open for discussion: cluster the columns?
    clust_cols <- hclust(dist(t(scale(weights_matrix)), method = "euclidean"),
                         method = "complete")
    
    # produce main heatmap plot
    hm <- pheatmap(mat = weights_matrix,
                   main = "",
                   color = colorRampPalette(c("white", "navy"))(50),
                   cutree_rows = 4,
                   cluster_rows = clust_rows,
                   cluster_cols = FALSE,
                   angle_col = 0,
                   treeheight_col = 20,
                   treeheight_row = 30,
                   legend = TRUE,
                   legend_breaks = c(0.01, 0.2, 0.4, 0.6, 0.8, 0.98, 0.99),
                   legend_labels = c(0,    0.2, 0.4, 0.6, 0.8, 1,    "AIC Weight\n\n"),
                   legend.position = "bottom")
    
    tiff(filename = file.path("./plots", exp_set, dosage, "AIC_heatmap.tiff"),
         height = 8, width = 6, units = 'in', res = 300, compression = "lzw")
    print(hm)
    dev.off()
    
    #########################
    # AIC Annotated Heatmap #
    #########################
    
    # produce an annotated heatmap
    hm_annotate <- data.frame(label = c("Log-Linear & Emax",
                                        "Mixed",
                                        "Quadratic & Emax",
                                        "Emax only"),
                              x = rep(-0.05, 4),
                              y = c(0.74, 0.415, 0.26, 0.115),
                              color = c("#000000", "#D55E00", "#000000", "#D55E00"))
    
    hm2 <- ggplotify::as.ggplot(hm$gtable) +
      geom_label(data = hm_annotate, aes(x = x, y = y, label = label, color = color),
                 hjust = "right", show.legend = FALSE) +
      scale_color_manual(values = c("#000000", "#D55E00")) +
      scale_x_continuous(limits = c(-0.3, 1)) +
      geom_errorbar(aes(x = -0.03, ymin = 0.52, ymax = 0.96), width = 0.02, lwd = 2, color = "#000000") +
      geom_errorbar(aes(x = -0.03, ymin = 0.32, ymax = 0.51), width = 0.02, lwd = 2, color = "#D55E00") +
      geom_errorbar(aes(x = -0.03, ymin = 0.21, ymax = 0.31), width = 0.02, lwd = 2, color = "#000000") +
      geom_errorbar(aes(x = -0.03, ymin = 0.03, ymax = 0.20), width = 0.02, lwd = 2, color = "#D55E00")
    
    tiff(filename = file.path("./plots", exp_set, dosage, "AIC_heatmap_ann.tiff"),
         height = 8, width = 8, units = 'in', res = 300, compression = "lzw")
    print(hm2)
    dev.off()
    
    ###################################
    # AIC Stacked, Clustered Bar Plot #
    ###################################
    
    # get dendrogram from study clusters
    dendro <- as.dendrogram(clust_rows)
    
    tiff(filename = file.path("./plots", exp_set, dosage, "barplot_weights_AIC_stack_clust.tiff"),
         height = 6, width = 8, units = 'in', res = 300, compression = "lzw")
    print(ggplot(data = weights_AIC %>%
                   # reverse order of dendrogram (to agree with convention)
                   mutate(trial_id = factor(trial_id, levels = rev(labels(dendro)))) %>%
                   # reverse order of columns
                   mutate(Model = factor(Model, levels = rev(levels(Model)))),
                 aes(fill = Model, x = AIC_weight, y = trial_id)) +
            theme_bw() +
            geom_bar(stat = "identity") +
            theme(axis.text.y = element_text(size = rel(0.9), angle = 0)) +
            # reverse colors to match reverse order of columns
            scale_fill_discrete(type = my_cols, limits = rev) +
            scale_color_discrete(type = my_cols) +
            scale_x_continuous(expand = c(0, 0), limits = c(-10e-10, 1+10e-10),
                               breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
            # reverse order of dendrogram
            ggh4x::scale_y_dendrogram(hclust = dendextend::rotate(clust_rows,
                                                                  order = rev(labels(dendro)),
                                                                  decreasing = TRUE)) +
            labs(x = "AIC Weight", y = "Trial ID") +
            geom_vline(xintercept = seq(0, 1, 0.2),
                       alpha = 0.1))
    dev.off()
    
    ################
    # AIC dot plot #
    ################
    
    tiff(filename = file.path("./plots", exp_set, dosage, "cleveland_dot_weights_AIC.tiff"),
         height = 6, width = 6, units = 'in', res = 300, compression = "lzw")
    set.seed(1)
    print(ggplot(weights_AIC, aes(x = AIC_weight, y = trial_id, color = Model)) +
            theme_bw() +
            geom_jitter(width = 0, height = 0, size = 4, alpha = 0.6) +
            scale_fill_discrete(type = my_cols) +
            scale_color_discrete(type = my_cols) +
            theme(axis.text.y = element_text(size = rel(0.9), angle = 0)) +
            labs(x = "AIC Weight", y = "Trial ID",
                 color = "Model",
                 caption = "Y-axis sorted by descending Emax weight"))
    dev.off()
    
    ####################################################
    # Max and nearly-min weights for each model family #
    ####################################################
    
    weights_AIC %>%
      group_by(Model) %>%
      select(-sum_sampsize) %>%
      filter(AIC_weight %in% max(AIC_weight))
    
    weights_AIC %>%
      group_by(Model) %>%
      select(-sum_sampsize) %>%
      filter(AIC_weight %in% c(min(AIC_weight), max(AIC_weight)) | AIC_weight < 0.15) %>%
      arrange(Model, AIC_weight)
    
    
    #########################
    # Calculate BIC weights #
    #########################
    
    # BIC weights
    weights_BIC <- data.frame(matrix(nrow = nrow(results_adj), ncol = 7))
    colnames(weights_BIC) <- c("trial_id", "sum_sampsize",
                               "eMax_weight", "linLog_weight",
                               "quad_weight", "lin_weight", 
                               "exp_weight")
    for (i in 1:nrow(results_adj)) {
      weights_BIC[i, 1:2] <- results_adj[i, c("trial_id", "sum_sampsize")]
      weights_BIC[i, 3:7] <- compute_weights(results_adj[i, c("eMax_BIC", "linLog_BIC",
                                                              "quad_BIC", "lin_BIC", 
                                                              "exp_BIC")])
    }
    
    # make numerical weight matrix
    weights_matrix <- weights_BIC[, 3:7]
    rownames(weights_matrix) <- weights_BIC$trial_id
    colnames(weights_matrix) <- c("Emax", "Log-Linear",
                                  "Quadratic", "Linear", 
                                  "Exponential")
    
    # wrangle weights into long format
    weights_BIC <- weights_BIC %>%
      arrange(desc(eMax_weight)) %>%
      mutate(trial_id = fct_reorder(trial_id, eMax_weight)) %>%
      pivot_longer(contains("_weight"), names_to = "Model", values_to = "BIC_weight") %>%
      mutate(Model = as.factor(as.character(Model))) %>%
      mutate(Model = recode_factor(Model,
                                   'eMax_weight'   = 'Emax',
                                   'linLog_weight' = 'Log-linear',
                                   'quad_weight'   = 'Quadratic',
                                   'lin_weight'    = 'Linear',
                                   'exp_weight'    = 'Exponential'))
    
    ##################
    # BIC Rank Table #
    ##################
    
    tab <- results_adj %>%
      as_tibble %>%
      select(eMax_BIC, linLog_BIC, lin_BIC, quad_BIC, exp_BIC) %>%
      as.matrix %>%
      rowRanks(x = ., ties.method = "max") %>%
      as_tibble %>%
      pivot_longer(everything(), names_to = "Model", values_to = "Rank") %>%
      mutate(Model = recode_factor(Model,
                                   'eMax_BIC' = 'Emax',
                                   'linLog_BIC' = 'Log-linear',
                                   'lin_BIC' = 'Linear',
                                   'quad_BIC' = 'Quadratic',
                                   'exp_BIC' = "Exponential")) %>%
      dplyr::count(Model, Rank) %>%
      pivot_wider(names_from = Rank, names_prefix = "Rank ",
                  values_from = n, values_fill = 0) %>%
      mutate(
        across(starts_with("Rank"), ~ ifelse(.x == sort(.x, partial = 5)[5], 0, 1), .names = "freq5_{col}"),
        across(starts_with("Rank"), ~ ifelse(.x == sort(.x, partial = 4)[4], 0, 1), .names = "freq4_{col}"),
        across(starts_with("Rank"), ~ ifelse(.x == sort(.x, partial = 3)[3], 0, 1), .names = "freq3_{col}"),
        across(starts_with("Rank"), ~ ifelse(.x == sort(.x, partial = 2)[2], 0, 1), .names = "freq2_{col}"),
        across(starts_with("Rank"), ~ ifelse(.x == sort(.x, partial = 1)[1], 0, 1), .names = "freq1_{col}")
      ) %>%
      arrange(pick(starts_with("freq"))) %>%
      select(Model, starts_with("Rank")) %>%
      mutate(y = 1:n()) %>%
      pivot_longer(cols = starts_with("Rank"), names_to = "Rank", values_to = "value") %>%
      mutate(x = -as.numeric(as.character(gsub("Rank ", "", Rank)))) %>%
      mutate(Rank = fct_reorder(Rank, -x),
             Model = fct_reorder(Model, -y))
    
    tiff(filename = file.path("./plots", exp_set, dosage, "model_ranks_BIC.tiff"),
         height = 2, width = 4, units = 'in', res = 300, compression = "lzw")
    print(ggplot(tab, aes(x = Rank, y = Model)) +
            geom_tile(aes(fill = value)) +
            geom_text(aes(label = value)) +
            scale_fill_gradient2(low = "white", 
                                 mid = "white", 
                                 high = "#228B22", 
                                 midpoint = 0) + # determine the colour
            theme(panel.grid.major.x = element_blank(), #no gridlines
                  panel.grid.minor.x = element_blank(), 
                  panel.grid.major.y = element_blank(), 
                  panel.grid.minor.y = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.ticks.y = element_blank(),
                  panel.background=element_rect(fill = "white"),
                  axis.text.x = element_text(size = 8, face = "bold"),
                  axis.text.y = element_text(size = 8, face = "bold")) + 
            theme(legend.title = element_text(face = "bold", size = 14),
                  legend.position = "none") + 
            scale_x_discrete(name = "", position = "top") +
            scale_y_discrete(name = ""))
    dev.off()
    
    ################
    # BIC dot plot #
    ################
    
    tiff(filename = file.path("./plots", exp_set, dosage, "cleveland_dot_weights_BIC.tiff"),
         height = 6, width = 6, units = 'in', res = 300, compression = "lzw")
    set.seed(1)
    print(ggplot(weights_BIC, aes(x = BIC_weight, y = trial_id, color = Model)) +
            theme_bw() +
            geom_jitter(width = 0, height = 0, size = 4, alpha = 0.6) +
            scale_fill_discrete(type = my_cols) +
            scale_color_discrete(type = my_cols) +
            theme(axis.text.y = element_text(size = rel(0.9), angle = 0)) +
            labs(x = "BIC Weight", y = "Trial ID",
                 color = "Model",
                 caption = "Y-axis sorted by descending Emax weight"))
    dev.off()
    
    #########################
    # Calculate BIC weights #
    #########################
    
    # BIC weights
    weights_BIC <- data.frame(matrix(nrow = nrow(results_adj), ncol = 7))
    colnames(weights_BIC) <- c("trial_id", "sum_sampsize",
                               "eMax_weight", "linLog_weight",
                               "quad_weight", "lin_weight", 
                               "exp_weight")
    for (i in 1:nrow(results_adj)) {
      weights_BIC[i, 1:2] <- results_adj[i, c("trial_id", "sum_sampsize")]
      weights_BIC[i, 3:7] <- compute_weights(results_adj[i, c("eMax_BIC", "linLog_BIC",
                                                              "quad_BIC", "lin_BIC", 
                                                              "exp_BIC")])
    }
    
    # make numerical weight matrix
    weights_matrix <- weights_BIC[, 3:7]
    rownames(weights_matrix) <- weights_BIC$trial_id
    colnames(weights_matrix) <- c("Emax", "Log-Linear",
                                  "Quadratic", "Linear", 
                                  "Exponential")
    
    # wrangle weights into long format
    weights_BIC <- weights_BIC %>%
      arrange(desc(eMax_weight)) %>%
      mutate(trial_id = fct_reorder(trial_id, eMax_weight)) %>%
      pivot_longer(contains("_weight"), names_to = "Model", values_to = "BIC_weight") %>%
      mutate(Model = as.factor(as.character(Model))) %>%
      mutate(Model = recode_factor(Model,
                                   'eMax_weight'   = 'Emax',
                                   'linLog_weight' = 'Log-linear',
                                   'quad_weight'   = 'Quadratic',
                                   'lin_weight'    = 'Linear',
                                   'exp_weight'    = 'Exponential'))
    
    ##################
    # BIC Rank Table #
    ##################
    
    tab <- results_adj %>%
      as_tibble %>%
      select(eMax_BIC, linLog_BIC, lin_BIC, quad_BIC, exp_BIC) %>%
      as.matrix %>%
      rowRanks(x = ., ties.method = "max") %>%
      as_tibble %>%
      pivot_longer(everything(), names_to = "Model", values_to = "Rank") %>%
      mutate(Model = recode_factor(Model,
                                   'eMax_BIC' = 'Emax',
                                   'linLog_BIC' = 'Log-linear',
                                   'lin_BIC' = 'Linear',
                                   'quad_BIC' = 'Quadratic',
                                   'exp_BIC' = "Exponential")) %>%
      dplyr::count(Model, Rank) %>%
      pivot_wider(names_from = Rank, names_prefix = "Rank ",
                  values_from = n, values_fill = 0) %>%
      mutate(
        across(starts_with("Rank"), ~ ifelse(.x == sort(.x, partial = 5)[5], 0, 1), .names = "freq5_{col}"),
        across(starts_with("Rank"), ~ ifelse(.x == sort(.x, partial = 4)[4], 0, 1), .names = "freq4_{col}"),
        across(starts_with("Rank"), ~ ifelse(.x == sort(.x, partial = 3)[3], 0, 1), .names = "freq3_{col}"),
        across(starts_with("Rank"), ~ ifelse(.x == sort(.x, partial = 2)[2], 0, 1), .names = "freq2_{col}"),
        across(starts_with("Rank"), ~ ifelse(.x == sort(.x, partial = 1)[1], 0, 1), .names = "freq1_{col}")
      ) %>%
      arrange(pick(starts_with("freq"))) %>%
      select(Model, starts_with("Rank")) %>%
      mutate(y = 1:n()) %>%
      pivot_longer(cols = starts_with("Rank"), names_to = "Rank", values_to = "value") %>%
      mutate(x = -as.numeric(as.character(gsub("Rank ", "", Rank)))) %>%
      mutate(Rank = fct_reorder(Rank, -x),
             Model = fct_reorder(Model, -y))
    
    tiff(filename = file.path("./plots", exp_set, dosage, "model_ranks_BIC.tiff"),
         height = 2, width = 4, units = 'in', res = 300, compression = "lzw")
    print(ggplot(tab, aes(x = Rank, y = Model)) +
            geom_tile(aes(fill = value)) +
            geom_text(aes(label = value)) +
            scale_fill_gradient2(low = "white", 
                                 mid = "white", 
                                 high = "#228B22", 
                                 midpoint = 0) + # determine the colour
            theme(panel.grid.major.x = element_blank(), #no gridlines
                  panel.grid.minor.x = element_blank(), 
                  panel.grid.major.y = element_blank(), 
                  panel.grid.minor.y = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.ticks.y = element_blank(),
                  panel.background=element_rect(fill = "white"),
                  axis.text.x = element_text(size = 8, face = "bold"),
                  axis.text.y = element_text(size = 8, face = "bold")) + 
            theme(legend.title = element_text(face = "bold", size = 14),
                  legend.position = "none") + 
            scale_x_discrete(name = "", position = "top") +
            scale_y_discrete(name = ""))
    dev.off()
    
    ################
    # BIC dot plot #
    ################
    
    tiff(filename = file.path("./plots", exp_set, dosage, "cleveland_dot_weights_BIC.tiff"),
         height = 6, width = 6, units = 'in', res = 300, compression = "lzw")
    set.seed(1)
    print(ggplot(weights_BIC, aes(x = BIC_weight, y = trial_id, color = Model)) +
            theme_bw() +
            geom_jitter(width = 0, height = 0, size = 4, alpha = 0.6) +
            scale_fill_discrete(type = my_cols) +
            scale_color_discrete(type = my_cols) +
            theme(axis.text.y = element_text(size = rel(0.9), angle = 0)) +
            labs(x = "BIC Weight", y = "Trial ID",
                 color = "Model",
                 caption = "Y-axis sorted by descending Emax weight"))
    dev.off()
    
    ############################################
    # Calculate AIC weights using Sigmoid-Emax #
    ############################################
    
    # AIC weights with Sigmoid-Emax
    weights_AIC <- data.frame(matrix(nrow = nrow(results_adj), ncol = 7))
    colnames(weights_AIC) <- c("trial_id", "sum_sampsize",
                               "sigEmax_weight", "linLog_weight",
                               "quad_weight", "lin_weight", 
                               "exp_weight")
    for (i in 1:nrow(results_adj)) {
      weights_AIC[i, 1:2] <- results_adj[i, c("trial_id", "sum_sampsize")]
      weights_AIC[i, 3:7] <- compute_weights(results_adj[i, c("sigEmax_AIC", "linLog_AIC",
                                                              "quad_AIC", "lin_AIC", 
                                                              "exp_AIC")])
    }
    
    # make numerical weight matrix
    weights_matrix <- weights_AIC[, 3:7]
    rownames(weights_matrix) <- weights_AIC$trial_id
    colnames(weights_matrix) <- c("Sigmoid-Emax", "Log-Linear",
                                  "Quadratic", "Linear", 
                                  "Exponential")
    
    # wrangle weights into long format
    weights_AIC <- weights_AIC %>%
      arrange(desc(sigEmax_weight)) %>%
      mutate(trial_id = fct_reorder(trial_id, sigEmax_weight)) %>%
      pivot_longer(contains("_weight"), names_to = "Model", values_to = "AIC_weight") %>%
      mutate(Model = as.factor(as.character(Model))) %>%
      mutate(Model = recode_factor(Model,
                                   'sigEmax_weight'   = 'Sigmoid-Emax',
                                   'linLog_weight' = 'Log-linear',
                                   'quad_weight'   = 'Quadratic',
                                   'lin_weight'    = 'Linear',
                                   'exp_weight'    = 'Exponential'))
    
    #####################################
    # AIC Rank Table using Sigmoid-Emax #
    #####################################
    
    tab <- results_adj %>%
      as_tibble %>%
      select(sigEmax_AIC, linLog_AIC, lin_AIC, quad_AIC, exp_AIC) %>%
      as.matrix %>%
      rowRanks(x = ., ties.method = "max") %>%
      as_tibble %>%
      pivot_longer(everything(), names_to = "Model", values_to = "Rank") %>%
      mutate(Model = recode_factor(Model,
                                   'sigEmax_AIC' = 'Sigmoid-Emax',
                                   'linLog_AIC' = 'Log-linear',
                                   'lin_AIC' = 'Linear',
                                   'quad_AIC' = 'Quadratic',
                                   'exp_AIC' = "Exponential")) %>%
      dplyr::count(Model, Rank) %>%
      pivot_wider(names_from = Rank, names_prefix = "Rank ",
                  values_from = n, values_fill = 0) %>%
      mutate(
        across(starts_with("Rank"), ~ ifelse(.x == sort(.x, partial = 5)[5], 0, 1), .names = "freq5_{col}"),
        across(starts_with("Rank"), ~ ifelse(.x == sort(.x, partial = 4)[4], 0, 1), .names = "freq4_{col}"),
        across(starts_with("Rank"), ~ ifelse(.x == sort(.x, partial = 3)[3], 0, 1), .names = "freq3_{col}"),
        across(starts_with("Rank"), ~ ifelse(.x == sort(.x, partial = 2)[2], 0, 1), .names = "freq2_{col}"),
        across(starts_with("Rank"), ~ ifelse(.x == sort(.x, partial = 1)[1], 0, 1), .names = "freq1_{col}")
      ) %>%
      arrange(pick(starts_with("freq"))) %>%
      select(Model, starts_with("Rank")) %>%
      mutate(y = 1:n()) %>%
      pivot_longer(cols = starts_with("Rank"), names_to = "Rank", values_to = "value") %>%
      mutate(x = -as.numeric(as.character(gsub("Rank ", "", Rank)))) %>%
      mutate(Rank = fct_reorder(Rank, -x),
             Model = fct_reorder(Model, -y))
    
    tiff(filename = file.path("./plots", exp_set, dosage, "model_ranks_AIC_se.tiff"),
         height = 2, width = 4, units = 'in', res = 300, compression = "lzw")
    print(ggplot(tab, aes(x = Rank, y = Model)) +
            geom_tile(aes(fill = value)) +
            geom_text(aes(label = value)) +
            scale_fill_gradient2(low = "white", 
                                 mid = "white", 
                                 high = "#228B22", 
                                 midpoint = 0) +
            theme(panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank(), 
                  panel.grid.major.y = element_blank(), 
                  panel.grid.minor.y = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.ticks.y = element_blank(),
                  panel.background=element_rect(fill = "white"),
                  axis.text.x = element_text(size = 8, face = "bold"),
                  axis.text.y = element_text(size = 8, face = "bold")) + 
            theme(legend.title = element_text(face = "bold", size = 14),
                  legend.position = "none") + 
            scale_x_discrete(name = "", position = "top") +
            scale_y_discrete(name = ""))
    dev.off()
    
    ###################
    # Waterfall plots #
    ###################
    
    # Sigmoid-Emax versus Hyperbolic-Emax
    results_adj <- results_adj %>%
      arrange(del_AIC_sigEmax) %>%
      mutate(order = 1:n() - 1) %>%
      mutate(sigEmax_better = factor(sigEmax_better,
                                     levels = c("Sigmoid-Emax", "Hyperbolic-Emax")))
    
    tiff(filename = file.path("./plots", exp_set, dosage, "waterfall_del_AIC_sigemax.tiff"),
         height = 3, width = 4, units = 'in', res = 300, compression = "lzw")
    print(ggplot(data = results_adj, aes(x = -order, y = del_AIC_sigEmax,
                                         color = sigEmax_better)) +
            geom_bar(stat = "identity", width = 1, orientation = "x", fill = "#FFFFFF",
                     linewidth = 0.5) +
            geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
            geom_hline(yintercept = c(-6, -2.2, 2.2, 6), lty = 2, color = "black",
                       linewidth = 0.5, alpha = 0.4) +
            geom_vline(xintercept = median(-results_adj$order), lty = 2, color = "black",
                       linewidth = 0.5, alpha = 0.4) +
            scale_color_discrete(name = "Winner: ", type = c("#E69F00", "#000000")) +
            scale_y_continuous(breaks = c(-6, -2.2, 0, 2.2, 6),
                               limits = c(-max(c(abs(results_adj$del_AIC_sigEmax), 6)),
                                          max(c(abs(results_adj$del_AIC_sigEmax), 6)))) +
            labs(x = NULL,
                 y = expression(Delta*" AIC")) +
            theme_classic() +
            theme(axis.line.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.y = element_text(face = "bold", angle = 90),
                  legend.position = "right"))
    dev.off()
    
    # Log-Linear versus Emax
    results_adj <- results_adj %>%
      arrange(del_AIC_linLog) %>%
      mutate(order = 1:n() - 1)
    
    tiff(filename = file.path("./plots", exp_set, dosage, "waterfall_del_AIC_linLog.tiff"),
         height = 3, width = 4, units = 'in', res = 300, compression = "lzw")
    print(ggplot(data = results_adj, aes(x = -order, y = del_AIC_linLog,
                                         color = linLog_better)) +
            geom_bar(stat = "identity", width = 1, orientation = "x", fill = "#FFFFFF",
                     linewidth = 0.5) +
            geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
            geom_hline(yintercept = c(-6, -2.2, 2.2, 6), lty = 2, color = "black",
                       linewidth = 0.5, alpha = 0.4) +
            geom_vline(xintercept = median(-results_adj$order), lty = 2, color = "black",
                       linewidth = 0.5, alpha = 0.4) +
            scale_color_discrete(name = "Winner: ", type = c("#000000", "#56B4E9")) +
            scale_y_continuous(breaks = c(-6, -2.2, 0, 2.2, 6),
                               limits = c(-max(c(abs(results_adj$del_AIC_linLog), 6)),
                                          max(c(abs(results_adj$del_AIC_linLog), 6)))) +
            labs(x = NULL,
                 y = expression(Delta*" AIC")) +
            theme_classic() +
            theme(axis.line.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.y = element_text(face = "bold", angle = 90),
                  legend.position = "right"))
    dev.off()
    
    # Linear versus Emax
    results_adj <- results_adj %>%
      arrange(del_AIC_lin) %>%
      mutate(order = 1:n() - 1)
    
    tiff(filename = file.path("./plots", exp_set, dosage, "waterfall_del_AIC_lin.tiff"),
         height = 2, width = 3, units = 'in', res = 300, compression = "lzw")
    print(ggplot(data = results_adj, aes(x = -order, y = del_AIC_lin,
                                         fill = lin_better, color = lin_better)) +
            geom_bar(stat = "identity", width = 1, orientation = "x") +
            geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
            geom_hline(yintercept = c(-6, -2.2, 2.2, 6), lty = 2, color = "black",
                       linewidth = 0.5, alpha = 0.4) +
            geom_vline(xintercept = median(-results_adj$order), lty = 2, color = "black",
                       linewidth = 0.5, alpha = 0.4) +
            scale_fill_discrete(name = "Winner", type = c("#000000", "#F0E442")) +
            scale_color_discrete(guide = "none", type = c("#000000", "#F0E442")) +
            scale_y_continuous(breaks = c(-6, -2.2, 0, 2.2, 6),
                               limits = c(-max(abs(results_adj$del_AIC_lin)),
                                          max(abs(results_adj$del_AIC_lin)))) +
            labs(x = NULL,
                 y = expression(Delta*" AIC")) +
            theme_classic() +
            theme(axis.line.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.y = element_text(face = "bold", angle = 90),
                  legend.position = "bottom",
                  #legend.title = element_blank(),
                  text = element_text(size = 15)))
    dev.off()
    
    # Quadratic versus Emax
    results_adj <- results_adj %>%
      arrange(del_AIC_quad) %>%
      mutate(order = 1:n() - 1)
    
    tiff(filename = file.path("./plots", exp_set, dosage, "waterfall_del_AIC_quad.tiff"),
         height = 2, width = 3, units = 'in', res = 300, compression = "lzw")
    print(ggplot(data = results_adj, aes(x = -order, y = del_AIC_quad,
                                         fill = quad_better, color = quad_better)) +
            geom_bar(stat = "identity", width = 1, orientation = "x") +
            geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
            geom_hline(yintercept = c(-6, -2.2, 2.2, 6), lty = 2, color = "black",
                       linewidth = 0.5, alpha = 0.4) +
            geom_vline(xintercept = median(-results_adj$order), lty = 2, color = "black",
                       linewidth = 0.5, alpha = 0.4) +
            scale_fill_discrete(name = "Winner", type = c("#000000", "#CC79A7")) +
            scale_color_discrete(guide = "none", type = c("#000000", "#CC79A7")) +
            scale_y_continuous(breaks = c(-6, -2.2, 0, 2.2, 6),
                               limits = c(-max(abs(results_adj$del_AIC_quad)),
                                          max(abs(results_adj$del_AIC_quad)))) +
            labs(x = NULL,
                 y = expression(Delta*" AIC")) +
            theme_classic() +
            theme(axis.line.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.y = element_text(face = "bold", angle = 90),
                  legend.position = "bottom",
                  #legend.title = element_blank(),
                  text = element_text(size = 15)))
    dev.off()
    
    # Exponential versus Emax
    results_adj <- results_adj %>%
      arrange(del_AIC_exp) %>%
      mutate(order = 1:n() - 1)
    
    tiff(filename = file.path("./plots", exp_set, dosage, "waterfall_del_AIC_exp.tiff"),
         height = 2, width = 3, units = 'in', res = 300, compression = "lzw")
    print(ggplot(data = results_adj, aes(x = -order, y = del_AIC_exp,
                                         fill = exp_better, color = exp_better)) +
            geom_bar(stat = "identity", width = 1, orientation = "x") +
            geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
            geom_hline(yintercept = c(-6, -2.2, 2.2, 6), lty = 2, color = "black",
                       linewidth = 0.5, alpha = 0.4) +
            geom_vline(xintercept = median(-results_adj$order), lty = 2, color = "black",
                       linewidth = 0.5, alpha = 0.4) +
            scale_fill_discrete(name = "Winner", type = c("#000000", "#009E73")) +
            scale_color_discrete(guide = "none", type = c("#000000", "#009E73")) +
            scale_y_continuous(breaks = c(-6, -2.2, 0, 2.2, 6),
                               limits = c(-max(abs(results_adj$del_AIC_exp)),
                                          max(abs(results_adj$del_AIC_exp)))) +
            labs(x = NULL,
                 y = expression(Delta*" AIC")) +
            theme_classic() +
            theme(axis.line.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.y = element_text(face = "bold", angle = 90),
                  legend.position = "bottom",
                  #legend.title = element_blank(),
                  text = element_text(size = 15)))
    dev.off()
    
    #################
    # Special cases #
    #################
    
    # Sigmoid-Emax beats Emax
    results_adj %>%
      arrange(desc(del_AIC_sigEmax)) %>%
      select(trial_id, del_AIC_sigEmax) %>%
      head(3)
    
    # Show id1007_1 and id1015_1 have identical data
    data(metaData)
    data_id1007_1 <- metaData %>%
      filter(taid == 1007 & protid == 1)
    data_id1015_1 <- metaData %>%
      filter(taid == 1015 & protid == 1)
    
    # identify columns that differ between id1007_1 and id1015_1
    differ <- names(which(!sapply(colnames(metaData),
                                  function(col){identical(data_id1007_1[, col],
                                                          data_id1015_1[, col])})))
    data_id1007_1 %>%
      select(all_of(differ))
    
    data_id1015_1 %>%
      select(all_of(differ))
  }
}
