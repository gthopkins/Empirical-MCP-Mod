# Author: Grant Hopkins
# Date: August 22, 2023
#
# Description: the script models mean of the response for continuous traits via
#              five parametric MCP-Mod models, in addition to an intercept-only
#              model, currently using (extra) weighted least squares linear or 
#              non-linear regression. Then the script calculates AIC and BIC for
#              each model.
#
# Input: a single tibble from md_list in the wrangle_data output
#
# Output: a list of lists for each model, where each model yields a list of
#         1) the model object, 2) the AIC, and 3) the BIC

# now constructed fitted models to the data
make_fits_cont <- function(data, default_exp = TRUE) {
  if(require(DoseFinding)){"You must library the DoseFinding package."}
  
  # extract essential information from data
  dose = data$dose
  sampsize = data$sampsize
  se = data$se
  rslt = data$rslt
  
  ##########
  # linLog #
  ##########
  
  linLog_mod_naive <- fitMod(dose = dose, resp = rslt, data = data,
                             model = "linlog", type = "general",
                             S = diag(data$se^2))
  
  linLog_naive <- coef(linLog_mod_naive)
  
  f_opt_linLog <- function(par) {
    e0 = par["e0"]
    delta = par["delta"]
    m_theta <- linlog(dose = dose, e0 = e0, delta = delta, off = 0.01 * max(dose))
    return(sum(sampsize * log((sampsize - 1) * se^2 + (rslt - m_theta)^2)))
  }
  
  linLog_opt <- optim(par = linLog_naive,
                      fn = f_opt_linLog,
                      method = "Nelder-Mead")
  
  linLog_MLE <- linLog_opt$par
  
  linLog_pred_MLE <- linlog(dose = data$dose,
                            e0 = linLog_MLE["e0"],
                            delta = linLog_MLE["delta"],
                            off = 0.01 * max(data$dose))
  
  linLog_maxLL <- -0.5 * (sum(sampsize * log((sampsize - 1) * se^2 + (rslt - linLog_pred_MLE)^2)) +
    sum(sampsize) * (log(2 * pi) + 1))
  
  linLog_AIC <- -2 * linLog_maxLL + length(linLog_MLE) * 2
  
  linLog_BIC <- -2 * linLog_maxLL + length(linLog_MLE) * log(sum(sampsize))
  
  #############
  # quadratic #
  #############
  
  quad_mod_naive <- fitMod(dose = dose, resp = rslt, data = data,
                           model = "quadratic", type = "general",
                           S = diag(data$se^2))
  
  quad_naive <- coef(quad_mod_naive)
  
  f_opt_quad <- function(par) {
    e0 = par["e0"]
    b1 = par["b1"]
    b2 = par["b2"]
    m_theta <- quadratic(dose = dose, e0 = e0, b1 = b1, b2 = b2)
    return(sum(sampsize * log((sampsize - 1) * se^2 + (rslt - m_theta)^2)))
  }
  
  quad_opt <- optim(par = quad_naive,
                    fn = f_opt_quad,
                    method = "Nelder-Mead")
  
  quad_MLE <- quad_opt$par
  
  quad_pred_MLE <- quadratic(dose = data$dose,
                             e0 = quad_MLE["e0"],
                             b1 = quad_MLE["b1"],
                             b2 = quad_MLE["b2"])
  
  quad_maxLL <- -0.5 * (sum(sampsize * log((sampsize - 1) * se^2 + (rslt - quad_pred_MLE)^2)) +
    sum(sampsize) * (log(2 * pi) + 1))
  
  quad_AIC <- -2 * quad_maxLL + length(quad_MLE) * 2
  
  quad_BIC <- -2 * quad_maxLL + length(quad_MLE) * log(sum(sampsize))
  
  ##########
  # linear #
  ##########
  
  lin_mod_naive <- fitMod(dose = dose, resp = rslt, data = data,
                          model = "linear", type = "general",
                          S = diag(data$se^2))
  
  lin_naive <- coef(lin_mod_naive)
  
  f_opt_lin <- function(par) {
    e0 = par["e0"]
    delta = par["delta"]
    m_theta <- linear(dose = dose, e0 = e0, delta = delta)
    return(sum(sampsize * log((sampsize - 1) * se^2 + (rslt - m_theta)^2)))
  }
  
  lin_opt <- optim(par = lin_naive,
                   fn = f_opt_lin,
                   method = "Nelder-Mead")
  
  lin_MLE <- lin_opt$par
  
  lin_pred_MLE <- linear(dose = data$dose,
                         e0 = lin_MLE["e0"],
                         delta = lin_MLE["delta"])
  
  lin_maxLL <- -0.5 * (sum(sampsize * log((sampsize - 1) * se^2 + (rslt - lin_pred_MLE)^2)) +
    sum(sampsize) * (log(2 * pi) + 1))
  
  lin_AIC <- -2 * lin_maxLL + length(lin_MLE) * 2
  
  lin_BIC <- -2 * lin_maxLL + length(lin_MLE) * log(sum(sampsize))
  
  #############
  # intercept #
  #############
  
  intercept_mod_naive <- lm(formula = rslt ~ 1,
                            data = data,
                            weights = data$se^2)
  
  intercept_naive <- coef(intercept_mod_naive)
  
  f_opt_intercept <- function(par) {
    e0 = par
    m_theta <- linear(dose = dose, e0 = e0, delta = 0)
    return(sum(sampsize * log((sampsize - 1) * se^2 + (rslt - m_theta)^2)))
  }
  
  intercept_opt <- optim(par = intercept_naive,
                         fn = f_opt_intercept,
                         method = "Brent",
                         lower = min(rslt),
                         upper = max(rslt))
  
  intercept_MLE <- setNames(intercept_opt$par, "e0")
  
  intercept_pred_MLE <- linear(dose = data$dose,
                               e0 = intercept_MLE["e0"],
                               delta = 0)
  
  intercept_maxLL <- -0.5 * (sum(sampsize * log((sampsize - 1) * se^2 + (rslt - intercept_pred_MLE)^2)) +
    sum(sampsize) * (log(2 * pi) + 1))
  
  intercept_AIC <- -2 * intercept_maxLL + length(intercept_MLE) * 2
  
  intercept_BIC <- -2 * intercept_maxLL + length(intercept_MLE) * log(sum(sampsize))
  
  ########
  # eMax #
  ########
  
  eMax_mod_naive <- fitMod(dose = dose, resp = rslt, data = data,
                           model = "emax", type = "general",
                           S = diag(data$se^2),
                           bnds = defBnds(max(data$dose))$emax)
  
  eMax_naive <- coef(eMax_mod_naive)
  r1 <- diff(range(defBnds(max(dose))$emax))
  if (eMax_naive["ed50"] >= defBnds(max(dose))$emax[2] - r1 / 1e3) {
    eMax_naive["ed50"] <- defBnds(max(dose))$emax[2] - r1 / 1e3} 
  if (eMax_naive["ed50"] <= defBnds(max(dose))$emax[1] + r1 / 1e3) {
    eMax_naive["ed50"] <- defBnds(max(dose))$emax[1] + r1 / 1e3}
  
  # eMax DoseFinding constraints
  cmatrix <- matrix(c(0, 0, +1,
                      0, 0, -1),
                    nrow = 2, byrow = TRUE)
  
  cci <- matrix(c(+defBnds(max(data$dose))$emax[1],
                  -defBnds(max(data$dose))$emax[2]))
  
  f_opt_eMax <- function(par) {
    e0 = par["e0"]
    eMax = par["eMax"]
    ed50 = par["ed50"]
    m_theta <- emax(dose = dose, e0 = e0, eMax = eMax, ed50 = ed50)
    return(sum(sampsize * log((sampsize - 1) * se^2 + (rslt - m_theta)^2)))
  }
  
  eMax_opt <- constrOptim(theta = eMax_naive,
                          f = f_opt_eMax,
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
  
  eMax_pred_MLE <- emax(dose = data$dose,
                        e0 = eMax_MLE["e0"],
                        eMax = eMax_MLE["eMax"],
                        ed50 = eMax_MLE["ed50"])
  
  eMax_maxLL <- -0.5 * (sum(sampsize * log((sampsize - 1) * se^2 + (rslt - eMax_pred_MLE)^2)) +
    sum(sampsize) * (log(2 * pi) + 1))
  
  eMax_AIC <- -2 * eMax_maxLL + length(eMax_MLE) * 2
  
  eMax_BIC <- -2 * eMax_maxLL + length(eMax_MLE) * log(sum(sampsize))
  
  ###############
  # exponential #
  ###############
  
  if (default_exp) {
    exp_bounds <- defBnds(max(dose))$exponential
  } else {
    exp_bounds <- c(-2, 2) * max(data$dose)
  }
  
  # manually deal with Inf values in grid search that result from dividing by 0
  if (data$trial_id[1] %in% c("id2_1", "id3_1", "id4007_1")){
    exp_mod_naive <- fitMod(dose = dose, resp = rslt, data = data,
                            model = "exponential", type = "general",
                            S = diag(length(unique(data$dose)) * data$se^2),
                            bnds = exp_bounds,
                            control = list(gridSize = list(dim1 = 100))) # default is 30
  } else {
    exp_mod_naive <- fitMod(dose = dose, resp = rslt, data = data,
                            model = "exponential", type = "general",
                            S = diag(length(unique(data$dose)) * data$se^2),
                            bnds = exp_bounds)
  }
  
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
  
  f_opt_exp <- function(par) {
    e0 = par["e0"]
    e1 = par["e1"]
    delta = par["delta"]
    m_theta <- exponential(dose = dose, e0 = e0, e1 = e1, delta = delta)
    return(sum(sampsize * log((sampsize - 1) * se^2 + (rslt - m_theta)^2)))
  }
  
  exp_opt <- constrOptim(theta = exp_naive,
                         f = f_opt_exp,
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
  
  exp_pred_MLE <- exponential(dose = data$dose,
                              e0 = exp_MLE["e0"],
                              e1 = exp_MLE["e1"],
                              delta = exp_MLE["delta"])
  
  exp_maxLL <- -0.5 * (sum(sampsize * log((sampsize - 1) * se^2 + (rslt - exp_pred_MLE)^2)) +
    sum(sampsize) * (log(2 * pi) + 1))
  
  exp_AIC <- -2 * exp_maxLL + length(exp_MLE) * 2
  
  exp_BIC <- -2 * exp_maxLL + length(exp_MLE) * log(sum(sampsize))
  
  ###########
  # sigEmax #
  ###########
  
  sigEmax_mod_naive <- fitMod(dose = dose, resp = rslt, data = data,
                              model = "sigEmax", type = "general",
                              S = diag(data$se^2),
                              bnds = defBnds(max(data$dose))$sigEmax)
  
  sigEmax_naive <- coef(sigEmax_mod_naive)
  
  r3 <- diff(range(defBnds(max(dose))$sigEmax[1, 1:2]))
  if (sigEmax_naive["ed50"] >= defBnds(max(dose))$sigEmax[1, 2] - r3 / 1e3) {
    sigEmax_naive["ed50"] <- defBnds(max(dose))$sigEmax[1, 2] - r3 / 1e3} 
  if (sigEmax_naive["ed50"] <= defBnds(max(dose))$sigEmax[1, 1] + r3 / 1e3) {
    sigEmax_naive["ed50"] <- defBnds(max(dose))$sigEmax[1, 1] + r3 / 1e3}
  
  r4 <- diff(range(defBnds(max(dose))$sigEmax[2, 1:2]))
  if (sigEmax_naive["h"] >= defBnds(max(dose))$sigEmax[2, 2] - r4 / 1e3) {
    sigEmax_naive["h"] <- defBnds(max(dose))$sigEmax[2, 2] - r4 / 1e3}
  if (sigEmax_naive["h"] <= defBnds(max(dose))$sigEmax[2, 1] + r4 / 1e3) {
    sigEmax_naive["h"] <- defBnds(max(dose))$sigEmax[2, 1] + r4 / 1e3}
  
  # sigEmax DoseFinding constraints
  cmatrix <- matrix(c(0, 0, +1,  0,
                      0, 0, -1,  0,
                      0, 0,  0, +1,
                      0, 0,  0, -1),
                    nrow = 4, byrow = TRUE)
  
  cci <- matrix(c(+defBnds(max(dose))$sigEmax[1, 1],
                  -defBnds(max(dose))$sigEmax[1, 2],
                  +defBnds(max(dose))$sigEmax[2, 1],
                  -defBnds(max(dose))$sigEmax[2, 2]))
  
  f_opt_sigEmax <- function(par) {
    e0 = par["e0"]
    eMax = par["eMax"]
    ed50 = par["ed50"]
    h = par["h"]
    m_theta <- sigEmax(dose = dose, e0 = e0, eMax = eMax, ed50 = ed50, h = h)
    return(sum(sampsize * log((sampsize - 1) * se^2 + (rslt - m_theta)^2)))
  }
  
  sigEmax_opt <- constrOptim(theta = sigEmax_naive,
                             f = f_opt_sigEmax,
                             grad = NULL,
                             method = "Nelder-Mead",
                             ui = cmatrix,
                             ci = cci,
                             mu = 1e-6,
                             control = list(),
                             outer.iterations = 1e3,
                             outer.eps = 1e-05,
                             hessian = FALSE)
  
  sigEmax_MLE <- sigEmax_opt$par
  
  sigEmax_pred_MLE <- sigEmax(dose = data$dose,
                              e0 = sigEmax_MLE["e0"],
                              eMax = sigEmax_MLE["eMax"],
                              ed50 = sigEmax_MLE["ed50"],
                              h = sigEmax_MLE["h"])
  
  sigEmax_maxLL <- -0.5 * (sum(sampsize * log((sampsize - 1) * se^2 + (rslt - sigEmax_pred_MLE)^2)) +
    sum(sampsize) * (log(2 * pi) + 1))
  
  sigEmax_AIC <- -2 * sigEmax_maxLL + length(sigEmax_MLE) * 2
  
  sigEmax_BIC <- -2 * sigEmax_maxLL + length(sigEmax_MLE) * log(sum(sampsize))
  
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
