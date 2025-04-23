# Overview of Empirical MCP-Mod code

In order to improve readability of code, I wrote functions to perform sub-tasks, as opposed to making a single file that needs to accomplish every task. The R-Markdown file **clinDR_analysis.Rmd** is the main file used to gather the results. The .R files produce functions that perform the sub-tasks. The **anim_pres_plots.Rmd** file is not pertinent to the analysis, but I used it in the presentation to make animated plots and am including it due to expressed interest. Now I will proceed to describe the sub-tasks:

## 1. pool_sum_stat.R

The **pool_sum_stat** function pools summary statistics when two dosing regimens in the same trial correspond to the same Total Daily Dose. 

For example, consider trial id1008_1. After the data processing steps above, it turns out that the total daily dose of 62.5 was given by two separate regimens, but the regimen is not deemed substantially different by the variable "primregimen" in the clinDR package. Hence, we assume regimen is independent from the response and we calculate summary statistics as if the two regimens were originally summarized together as follows:

|dose | regimen | sampsize| rslt  | sd    | se    |
|:---:|:-------:|:-------:|:-----:|:-----:|:-----:|
|62.5 | QD      | 59      | 0.049 | 0.169 | 0.022 |
|62.5 | BID     | 57      | 0.065 | 0.174 | 0.023 |

and returning:

|dose | regimen | sampsize| rslt  | sd    | se    |
|:---:|:-------:|:-------:|:-----:|:-----:|:-----:|
|62.5 | *N/A*   | 116     | 0.057 | 0.171 | 0.016 | 

where:  
$\qquad \small \text{sampsize} := n_p = 59 + 57 = 116$  

$\qquad \small \text{rslt} := \hat{\mu}_p = 0.049 \times \dfrac{59}{n_p} + 0.065 \times \dfrac{57}{n_p} = 0.057$  

$\qquad \small \text{sd} := \hat{\sigma}_p = \sqrt{\dfrac{(59 - 1) * 0.169^2 + 59 * (0.049 - \hat{\mu}_p)^2 + (57 - 1) * 0.174^2 + 57 * (0.065 - \hat{\mu}_p)^2}{n_p - 1}} = 0.171$

$\qquad \small \text{se} := \text{SE}_{\hat{\mu}_p} = \hat{\sigma}_p / \sqrt{n_p} = 0.016$  

Under the assumption that regimen is independent of the response, then the combined summary statistics are better estimators for the underlying parameters. Under our normality assumption, we must use the combined standard deviation in order to determine the maximum liklieood estimators for $\theta$ in $m(\omega;\theta)$. This explains why we do not simply find the best fitting parameters using two different means and standard errors at the same dose level.

## 2. wrangle_data.R

The **wrangle_data** function takes in the *metaData* dataset in the clinDR package. It returns processed data based on our inclusion criteria in two forms: 1) a list of study-level datasets and 2) a master dataset with study-level IDs.

## 3. validate_data.R

The **validate_data** function takes in one study-level dataset. It performs basic checks to make sure it looks as expected and returns an error if there are unexpected issues.

## 4. make_fits.R

The **make_fits** function takes in one study-level dataset. It returns a list of fitted model objects for the MCP-Mod parametric models, in addition to calculations of deviance, AIC, and BIC.

## 5. make_plot.R

The **make_plot** function takes in one study-level dataset and one study-level list of model fits. It ouputs a $2 \times 3$ grid of the raw data plotted with 5 parametric fits and the intercept-only fit.

## 6. compute_weights.R

The **compute_weights** functions takes in one study-level vector of model information criteria. It outputs the model weights according to those information criteria in the same order.
