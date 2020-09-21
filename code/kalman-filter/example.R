#!/usr/bin/env R

options(repos = "https://cran.microsoft.com/snapshot/2020-02-29/")
library(pomp)

source("../model-setup/makepompmodel.R")
source("../model-setup/setparsvars.R")

pomp_data <- data.frame(time = c(1, 100), cases = NA, hosps = NA, deaths = NA)

n_knots <- 2
est_these_pars = c("log_sigma_dw", "min_frac_dead", "max_frac_dead", "log_half_dead",
                   "log_theta_cases", "log_theta_deaths")
est_these_inivals = c("E1_0", "Ia1_0", "Isu1_0", "Isd1_0")
# est_these_inivals = ""  # to not estimate any initial values
knot_coefs <-  paste0("b", 1:n_knots)
est_these_pars <- c(est_these_pars, knot_coefs)

# Set the parameter values and initial conditions
par_var_list <- setparsvars(est_these_pars = est_these_pars, 
                            est_these_inivals = est_these_inivals,
                            population = 1e6,
                            rnaught = 6)  # set R0 at beginning of epidemic

pomp_covar <- covariate_table(
  t = pomp_data$time,
  seas = bspline.basis(
    x=t,
    nbasis=n_knots,
    degree=0
  ),
  rel_beta_change = 1,
  trend_sim = 10,  # this is a placeholder only needed for simulation
  fit = 1,  # 1 = fitting; 0 = simulating
  times="t",
  order = "constant"
)

pomp_model <- makepompmodel(par_var_list = par_var_list, 
                            pomp_data = pomp_data, 
                            pomp_covar = pomp_covar,
                            n_knots = n_knots)

coef(pomp_model) <- par_var_list$allparvals
out <- pomp::simulate(pomp_model, times = 1:30)

sd <- as.data.frame(out)

