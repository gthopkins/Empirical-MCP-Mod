# Author: Grant Hopkins
# Date: August 24, 2023
#
# Description: the script models odds of the response for binary traits via five
#              parametric MCP-Mod models, in addition to an intercept-only
#              model, primarily using (extra) non-linear logistic regression.
#              Then the script calculates AIC for each model.
#
# Input: a single tibble from md_list in the wrangle_data output
#
# Output: a list of lists for each model, where each model yields a list of
#         1) the model object, 2) the AIC, and 3) the BIC

# make function to calculate count of samples with binary response
count_data <- function(data){
  data %>%
    mutate(n1 = round(rslt * sampsize),
           n0 = sampsize - n1) %>%
    rename(rslt_prop = rslt) %>%
    rename(sampsize_tot = sampsize) %>%
    pivot_longer(cols = c(n1, n0), names_to = "rslt", values_to = "sampsize") %>%
    mutate(rslt = as.numeric(as.character(gsub("n", "", rslt)))) %>%
    select(trial_id, dose, rslt, sampsize)
}

# define inverse of logit
inv_logit <- function(x){1 - 1 / (1 + exp(x))}

# now constructed fitted models to the data
make_fits_bin <- function(data, default_exp = TRUE) {
  if(require(DoseFinding)){"You must library the DoseFinding package."}
  
  # convert into convenient binary format
  count_data <- count_data(data)
  
  # extract essential information from data
  bdose = count_data$dose
  brslt = count_data$rslt
  bsampsize = count_data$sampsize
  
  ######################################################
  ######################################################
  ### logistic regression for non-linear mean models ###
  ######################################################
  ######################################################
  
  ##########
  # linLog #
  ##########
  
  h_shift <- 0.01 * max(bdose)
  
  linLog_mod <- glm(rslt ~ 1 + I(log(dose + h_shift)),
                    data = count_data,
                    family = binomial(link = "logit"),
                    weights = bsampsize)
  
  linLog_MLE <- setNames(coef(linLog_mod),
                         c("e0", "delta"))
  
  linLog_pred_MLE <- linlog(dose = bdose,
                            e0 = linLog_MLE["e0"],
                            delta = linLog_MLE["delta"],
                            off = h_shift)
  
  linLog_maxLL <- sum(brslt * bsampsize * log(inv_logit(linLog_pred_MLE)) +
                        (1 - brslt) * bsampsize * (log(1 - inv_logit(linLog_pred_MLE))))
  
  linLog_AIC <- -2 * linLog_maxLL + length(linLog_MLE) * 2
  
  linLog_BIC <- -2 * linLog_maxLL + length(linLog_MLE) * log(sum(bsampsize))
  
  #############
  # quadratic #
  #############
  
  quad_mod <- glm(rslt ~ 1 + dose + I(dose^2),
                  data = count_data,
                  family = binomial(link = "logit"),
                  weights = bsampsize)
  
  quad_MLE <- setNames(coef(quad_mod),
                       c("e0", "b1", "b2"))
  
  quad_pred_MLE <- quadratic(dose = bdose,
                             e0 = quad_MLE["e0"],
                             b1 = quad_MLE["b1"],
                             b2 = quad_MLE["b2"])
  
  quad_maxLL <- sum(brslt * bsampsize * log(inv_logit(quad_pred_MLE)) +
                      (1 - brslt) * bsampsize * (log(1 - inv_logit(quad_pred_MLE))))
  
  quad_AIC <- -2 * quad_maxLL + length(quad_MLE) * 2
  
  quad_BIC <- -2 * quad_maxLL + length(quad_MLE) * log(sum(bsampsize))
  
  ##########
  # linear #
  ##########
  
  lin_mod <- glm(rslt ~ 1 + dose,
                 data = count_data,
                 family = binomial(link = "logit"),
                 weights = bsampsize)
  
  lin_MLE <- setNames(coef(lin_mod),
                      c("e0", "delta"))
  
  lin_pred_MLE <- linear(dose = bdose,
                         e0 = lin_MLE["e0"],
                         delta = lin_MLE["delta"])
  
  lin_maxLL <- sum(brslt * bsampsize * log(inv_logit(lin_pred_MLE)) +
                     (1 - brslt) * bsampsize * (log(1 - inv_logit(lin_pred_MLE))))
  
  lin_AIC <- -2 * lin_maxLL + length(lin_MLE) * 2
  
  lin_BIC <- -2 * lin_maxLL + length(lin_MLE) * log(sum(bsampsize))
  
  #############
  # intercept #
  #############
  
  intercept_mod <- glm(rslt ~ 1,
                       data = count_data,
                       family = binomial(link = "logit"),
                       weights = bsampsize)
  
  intercept_MLE <- setNames(coef(intercept_mod),
                            c("e0"))
  
  intercept_pred_MLE <- linear(dose = bdose,
                               e0 = intercept_MLE["e0"],
                               delta = 0)
  
  intercept_maxLL <- sum(brslt * bsampsize * log(inv_logit(intercept_pred_MLE)) +
                           (1 - brslt) * bsampsize * (log(1 - inv_logit(intercept_pred_MLE))))
  
  intercept_AIC <- -2 * intercept_maxLL + length(intercept_MLE) * 2
  
  intercept_BIC <- -2 * intercept_maxLL + length(intercept_MLE) * log(sum(bsampsize))
  
  ######################################################
  ######################################################
  ### logistic regression for non-linear mean models ###
  ######################################################
  ######################################################
  
  # transform data for parameter initialization
  data <- data %>%
    # use asymptotic standard deviation of logit-transformed proportion
    mutate(sd = sqrt(1 / sampsize * rslt * (1 - rslt))) %>%
    # now calculate the standard error and logit-transformed proportion itself
    mutate(se = sd / sqrt(sampsize),
           rslt = log(rslt / (1 - rslt)))
  
  ########
  # eMax #
  ########
  
  eMax_mod_naive <- fitMod(dose = dose, resp = rslt, data = data,
                           model = "emax", type = "general",
                           S = diag(data$se^2),
                           bnds = defBnds(max(bdose))$emax)
  
  eMax_naive <- coef(eMax_mod_naive)
  r1 <- diff(range(defBnds(max(bdose))$emax))
  if (eMax_naive["ed50"] >= defBnds(max(bdose))$emax[2] - r1 / 1e3) {
    eMax_naive["ed50"] <- defBnds(max(bdose))$emax[2] - r1 / 1e3} 
  if (eMax_naive["ed50"] <= defBnds(max(bdose))$emax[1] + r1 / 1e3) {
    eMax_naive["ed50"] <- defBnds(max(bdose))$emax[1] + r1 / 1e3}
  
  # eMax DoseFinding constraints
  cmatrix <- matrix(c(0, 0, +1,
                      0, 0, -1),
                    nrow = 2, byrow = TRUE)
  
  cci <- matrix(c(+defBnds(max(bdose))$emax[1],
                  -defBnds(max(bdose))$emax[2]))
  
  eMax_nLL <- function(par){
    e0 = par["e0"]
    eMax = par["eMax"]
    ed50 = par["ed50"]
    pred <- emax(dose = bdose, e0 = e0, eMax = eMax, ed50 = ed50)
    LL <- sum(brslt * bsampsize * log(inv_logit(pred)) +
                (1 - brslt) * bsampsize * (log(1 - inv_logit(pred))))
    return(-LL)
  }
  
  eMax_opt <- constrOptim(theta = eMax_naive,
                          f = eMax_nLL,
                          grad = NULL,
                          method = "Nelder-Mead",
                          ui = cmatrix,
                          ci = cci,
                          mu = 1e-6,
                          control = list(),
                          outer.iterations = 1e3,
                          outer.eps = 1e-05,
                          hessian = FALSE)
  
  eMax_MLE <- eMax_opt$par
  
  eMax_pred_MLE <- emax(dose = bdose,
                        e0 = eMax_MLE["e0"],
                        eMax = eMax_MLE["eMax"],
                        ed50 = eMax_MLE["ed50"])
  
  eMax_maxLL <- sum(brslt * bsampsize * log(inv_logit(eMax_pred_MLE)) +
                      (1 - brslt) * bsampsize * (log(1 - inv_logit(eMax_pred_MLE))))
  
  eMax_AIC <- -2 * eMax_maxLL + length(eMax_MLE) * 2
  
  eMax_BIC <- -2 * eMax_maxLL + length(eMax_MLE) * log(sum(bsampsize))
  
  ###############
  # exponential #
  ###############
  
  if (default_exp) {
    exp_bounds <- defBnds(max(bdose))$exponential
  } else {
    exp_bounds <- c(-2, 2) * max(bdose)
  }
  
  exp_mod_naive <- fitMod(dose = dose, resp = rslt, data = data,
                          model = "exponential", type = "general",
                          S = diag(data$se^2),
                          bnds = exp_bounds)
  
  exp_naive <- coef(exp_mod_naive)
  r2 <- diff(range(exp_bounds))
  if (exp_naive["delta"] >= exp_bounds[2] - r2 / 1e3) {
    exp_naive["delta"] <- exp_bounds[2] - r2 / 1e3} 
  if (exp_naive["delta"] <= exp_bounds[1] + r2 / 1e3) {
    exp_naive["delta"] <- exp_bounds[1] + r2 / 1e3}
  
  # exp DoseFinding constraints
  cmatrix <- matrix(c(0, 0, +1,
                      0, 0, -1),
                    nrow = 2, byrow = TRUE)
  
  cci <- matrix(c(+exp_bounds[1],
                  -exp_bounds[2]))
  
  exp_nLL <- function (par){
    e0 = par["e0"]
    e1 = par["e1"]
    delta = par["delta"]
    pred <- exponential(dose = bdose, e0 = e0, e1 = e1, delta = delta)
    LL <- sum(brslt * bsampsize * log(inv_logit(pred)) +
                (1 - brslt) * bsampsize * (log(1 - inv_logit(pred))))
    return(-LL)
  }
  
  exp_opt <- constrOptim(theta = exp_naive,
                         f = exp_nLL,
                         grad = NULL,
                         method = "Nelder-Mead",
                         ui = cmatrix,
                         ci = cci,
                         mu = 1e-6,
                         control = list(),
                         outer.iterations = 1e3,
                         outer.eps = 1e-05,
                         hessian = FALSE)
  
  exp_MLE <- exp_opt$par
  
  exp_pred_MLE <- exponential(dose = bdose,
                              e0 = exp_MLE["e0"],
                              e1 = exp_MLE["e1"],
                              delta = exp_MLE["delta"])
  
  exp_maxLL <- sum(brslt * bsampsize * log(inv_logit(exp_pred_MLE)) +
                     (1 - brslt) * bsampsize * (log(1 - inv_logit(exp_pred_MLE))))
  
  exp_AIC <- -2 * exp_maxLL + length(exp_MLE) * 2
  
  exp_BIC <- -2 * exp_maxLL + length(exp_MLE) * log(sum(bsampsize))
  
  ###########
  # sigEmax #
  ###########
  
  sigEmax_mod_naive <- fitMod(dose = dose, resp = rslt, data = data,
                              model = "sigEmax", type = "general",
                              S = diag(data$se^2),
                              bnds = defBnds(max(bdose))$sigEmax)
  
  sigEmax_naive <- coef(sigEmax_mod_naive)
  
  r3 <- diff(range(defBnds(max(bdose))$sigEmax[1, 1:2]))
  if (sigEmax_naive["ed50"] >= defBnds(max(bdose))$sigEmax[1, 2] - r3 / 1e3) {
    sigEmax_naive["ed50"] <- defBnds(max(bdose))$sigEmax[1, 2] - r3 / 1e3} 
  if (sigEmax_naive["ed50"] <= defBnds(max(bdose))$sigEmax[1, 1] + r3 / 1e3) {
    sigEmax_naive["ed50"] <- defBnds(max(bdose))$sigEmax[1, 1] + r3 / 1e3}
  
  r4 <- diff(range(defBnds(max(bdose))$sigEmax[2, 1:2]))
  if (sigEmax_naive["h"] >= defBnds(max(bdose))$sigEmax[2, 2] - r4 / 1e3) {
    sigEmax_naive["h"] <- defBnds(max(bdose))$sigEmax[2, 2] - r4 / 1e3}
  if (sigEmax_naive["h"] <= defBnds(max(bdose))$sigEmax[2, 1] + r4 / 1e3) {
    sigEmax_naive["h"] <- defBnds(max(bdose))$sigEmax[2, 1] + r4 / 1e3}
  
  # sigEmax DoseFinding constraints
  cmatrix <- matrix(c(0, 0, +1,  0,
                      0, 0, -1,  0,
                      0, 0,  0, +1,
                      0, 0,  0, -1),
                    nrow = 4, byrow = TRUE)
  
  cci <- matrix(c(+defBnds(max(bdose))$sigEmax[1, 1],
                  -defBnds(max(bdose))$sigEmax[1, 2],
                  +defBnds(max(bdose))$sigEmax[2, 1],
                  -defBnds(max(bdose))$sigEmax[2, 2]))
  
  sigEmax_nLL <- function (par){
    e0 = par["e0"]
    eMax = par["eMax"]
    ed50 = par["ed50"]
    h = par["h"]
    pred <- sigEmax(dose = bdose, e0 = e0, eMax = eMax, ed50 = ed50, h = h)
    LL <- sum(brslt * bsampsize * log(inv_logit(pred)) +
                (1 - brslt) * bsampsize * (log(1 - inv_logit(pred))))
    return(-LL)
  }
  
  sigEmax_opt <- constrOptim(theta = sigEmax_naive,
                             f = sigEmax_nLL,
                             grad = NULL,
                             method = "Nelder-Mead",
                             ui = cmatrix,
                             ci = cci,
                             mu = 1e-6,
                             control = list(maxit = 1e3),
                             outer.iterations = 1e3,
                             outer.eps = 1e-05,
                             hessian = FALSE)
  
  sigEmax_MLE <- sigEmax_opt$par
  
  sigEmax_pred_MLE <- sigEmax(dose = bdose,
                              e0 = sigEmax_MLE["e0"],
                              eMax = sigEmax_MLE["eMax"],
                              ed50 = sigEmax_MLE["ed50"],
                              h = sigEmax_MLE["h"])
  
  sigEmax_maxLL <- sum(brslt * bsampsize * log(inv_logit(sigEmax_pred_MLE)) +
                         (1 - brslt) * bsampsize * (log(1 - inv_logit(sigEmax_pred_MLE))))
  
  sigEmax_AIC <- -2 * sigEmax_maxLL + length(sigEmax_MLE) * 2
  
  sigEmax_BIC <- -2 * sigEmax_maxLL + length(sigEmax_MLE) * log(sum(bsampsize))
  
  ###############################
  # Combine results into a list #
  ###############################
  
  fit_list <- list(eMax = list(MLE = eMax_MLE,
                               maxLL = eMax_maxLL,
                               AIC = eMax_AIC,
                               BIC = eMax_BIC),
                   linLog = list(MLE = linLog_MLE,
                                 maxLL = linLog_maxLL,
                                 AIC = linLog_AIC,
                                 BIC = linLog_BIC),
                   exp = list(MLE = exp_MLE,
                              maxLL = exp_maxLL,
                              AIC = exp_AIC,
                              BIC = exp_BIC),
                   quad = list(MLE = quad_MLE,
                               maxLL = quad_maxLL,
                               AIC = quad_AIC,
                               BIC = quad_BIC),
                   lin = list(MLE = lin_MLE,
                              maxLL = lin_maxLL,
                              AIC = lin_AIC,
                              BIC = lin_BIC),
                   intercept = list(MLE = intercept_MLE,
                                    maxLL = intercept_maxLL,
                                    AIC = intercept_AIC,
                                    BIC = intercept_BIC),
                   sigEmax = list(MLE = sigEmax_MLE,
                                  maxLL = sigEmax_maxLL,
                                  AIC = sigEmax_AIC,
                                  BIC = sigEmax_BIC))
  
  return(fit_list)
}
