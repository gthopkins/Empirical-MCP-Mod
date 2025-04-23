# Author: Grant Hopkins
# Date: August 24, 2023
#
# Description: the script produces a plot of the fitted values along with 
#              standard errors of the summary statistics for binary response
#
# Input: fitted models objects from make_fits_bin.R and its corresponding
#        tibble from md_list in the wrangle_data output
#
# Output: plot and data used to generate the plot

# produce plot
make_plot_bin <- function(fits, data) {
  
  ######################################################
  # the section below produces the predicted fits data #
  ######################################################
  
  # identify entire dose range
  dose_range <- seq(0, max(data$dose), max(data$dose) / 1e3)
  
  
  # make predictions using parametric least-squares estimates
  model_fits <- data.frame(dose = dose_range,
          eMax_fit = emax(dose  = dose_range,
                          e0    = fits$eMax$MLE["e0"],
                          eMax  = fits$eMax$MLE["eMax"],
                          ed50  = fits$eMax$MLE["ed50"]),
          linlog_fit = linlog(dose  = dose_range,
                              e0    = fits$linLog$MLE["e0"],
                              delta = fits$linLog$MLE["delta"],
                              off   = 0.01 * max(data$dose)),
          lin_fit = linear(dose  = dose_range,
                           e0    = fits$lin$MLE["e0"],
                           delta = fits$lin$MLE["delta"]),
          quad_fit = quadratic(dose  = dose_range,
                               e0    = fits$quad$MLE["e0"],
                               b1    = fits$quad$MLE["b1"],
                               b2    = fits$quad$MLE["b2"]),
          exp_fit = exponential(dose  = dose_range,
                                e0    = fits$exp$MLE["e0"],
                                e1    = fits$exp$MLE["e1"],
                                delta = fits$exp$MLE["delta"]),
          null_fit = linear(dose  = dose_range,
                            e0    = fits$intercept$MLE["e0"],
                            delta = 0)) %>%
    
    # now make into long format
    pivot_longer(cols = contains("_fit"),
                 names_to = "model", values_to = "fit") %>%
    
    # coerce the model identifier into a factor
    mutate(model = as.factor(as.character(model))) %>%
    
    # order the factor labels by preference
    mutate(model = recode_factor(model,
                                 'eMax_fit' = 'Emax',
                                 'linlog_fit' = 'Log-linear',
                                 'quad_fit' = 'Quadratic',
                                 'lin_fit' = 'Linear',
                                 'exp_fit' = "Exponential",
                                 'null_fit' = 'Intercept')) %>%
    
    # convert the fitted values to the probability scale
    mutate(fit = inv_logit(fit))
  
  ###############################################################
  # the section below finds a suitable scale to plot the y-axis #
  ###############################################################
  
  # # identify suitable upper- and lower-bounds for y-axis
  # y_axis_range <- c(0, 1)
  # y_axis_range <- y_axis_range + c(-0.12, 0.01) * diff(y_axis_range)
  
  # identify suitable upper- and lower-bounds for y-axis
  y_axis_range <- range(c(data$rslt - data$se,
                          data$rslt + data$se,
                          model_fits$fit))
  y_axis_range <- y_axis_range + c(-0.12, 0.01) * diff(y_axis_range)
  
  ####################################################
  # the section below formats the AIC and BIC labels #
  ####################################################
  
  # identify where to place the AIC and BIC labels on the x- and y-axis
  x_place <- max(data$dose) - 0.28 * diff(range(data$dose))
  y_place <- min(y_axis_range) + 0.06 * diff(y_axis_range)
  
  # write out AIC and BIC labels for each plot, with careful care for the order
  AIC_vec <- c(fits$eMax$AIC,
               fits$linLog$AIC,
               fits$quad$AIC,
               fits$lin$AIC,
               fits$exp$AIC,
               fits$intercept$AIC)
  
  BIC_vec <- c(fits$eMax$BIC,
               fits$linLog$BIC,
               fits$quad$BIC,
               fits$lin$BIC,
               fits$exp$BIC,
               fits$intercept$BIC)
  
  # determine largest number of (finite) integer digits to align numbers
  del_AIC <- AIC_vec[is.finite(AIC_vec)] - min(AIC_vec[is.finite(AIC_vec)])
  del_BIC <- BIC_vec[is.finite(BIC_vec)] - min(BIC_vec[is.finite(BIC_vec)])
  max_p10 <- c(del_AIC, del_BIC) %>% log10 %>% floor %>% max + 1
  
  # align the numbers at the decimal point
  AIC_align <- format(round(AIC_vec - min(AIC_vec), 1),
                      nsmall = 1, width = max_p10 + 2, justify = "right")
  BIC_align <- format(round(BIC_vec - min(BIC_vec), 1),
                      nsmall = 1, width = max_p10 + 2, justify = "right")
  
  # produce dataset to label AIC and BIC for each model
  ann_text <- data.frame(x_place = x_place,
                         y_place = y_place,
                         label = paste0(paste0("\u0394AIC: ", AIC_align, "\n"),
                                        paste0("\u0394BIC: ", BIC_align)),
                         model = c('Emax',
                                   'Log-linear',
                                   'Quadratic',
                                   'Linear',
                                   'Exponential',
                                   'Intercept')) %>%
    mutate(model = as.factor(as.character(model)))
  
  ###########################################
  # the section below makes the actual plot #
  ###########################################
  
  # identify plot-friendly colors
  my_cols <- c("#000000", "#56B4E9", "#CC79A7", "#F0E442",
               "#009E73", "#D55E00", "#0072B2", "#E69F00")
  
  plot <- ggplot(data = model_fits) +
    # format plot background
    theme_bw() +
    
    # plot the summary statistics from the clinDR package
    geom_point(data = data, aes(x = dose, y = rslt), cex = 4, color = "black") +
    geom_linerange(data = data, aes(x = dose, ymin = rslt - se, ymax = rslt + se),
                   lwd = 0.4, alpha = 0.6, color = "black") +
    
    # # plot the sample sizes at top and bottom
    # geom_point(data = summary_data,
    #            aes(x = dose, y = rslt, size = size), color = "steelblue") +
    
    # plot the fitted models
    geom_line(aes(x = dose, y = fit, color = model), lwd = 1) +
    
    # set colors for the models
    scale_fill_discrete(type = my_cols) +
    scale_color_discrete(type = my_cols) +
    
    # facet by the model type
    facet_wrap(vars(model), scales = "fixed") +
    
    # label the plot with summary information
    labs(x = "Dose", y = "Estimated Proportion",
         caption = "\u0394AIC/\u0394BIC is difference from minimum AIC/BIC") +
    
    # make small formatting changes
    theme(text = element_text(size = 15),
          legend.position = "none") +
    
    # set y-axis to suitable range
    scale_y_continuous(limits = y_axis_range) +

    # add AIC and BIC annotations to each panel
    geom_label(data = ann_text, aes(x = x_place, y = y_place, label = label),
               fill = "white", alpha = 0.6, label.size = NA, hjust = 0)
  
  return(list(plot = plot, plot_data = model_fits))
}
