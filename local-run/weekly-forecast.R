#! /usr/bin/env Rscript

# weekly-forecast.R
# This script is designed to simulate the model forward under different 
# scenarios to produce the data needed to generate forecasts in the standard
# covidhub weekly format.

# Start with a clean workspace to avoid downstream errors -----------------
rm(list = ls(all.names = TRUE))

# --------------------------------------------------
# Load necessary libraries
# --------------------------------------------------
# could be loaded here or inside functions
library('lubridate')
library('purrr')
library('pomp')
library('dplyr')
library('tidyr')
library('tibble')
library(readr)

# --------------------------------------------------
# Source all needed functions/scripts
# --------------------------------------------------

source("../code/forward-simulations/simulate_trajectories.R")
source("../code/forward-simulations/runscenarios.R") #runs forward simulations 
source("../code/model-setup/makepompmodel.R") #generates the pomp model
source("../code/model-setup/setparsvars_warm.R") #setting all parameters, specifying those that are  fitted


# Load the pomp information
locname <- Sys.getenv("locname") 
fdt <- Sys.getenv("fdt")

## load version of data used in model fitting
data_path <- file.path("archive", fdt, paste0(locname, ".csv"))

mdata <- read_csv(
  data_path,
  col_types = cols_only(
    location = col_character(),
    date = col_date(format = ""),
    variable = col_character(),
    mean_value = col_double()
  )
)

mdata2 <- mdata %>% 
  select(date, location, variable, mean_value) %>%
  filter(variable %in% c("actual_daily_cases", "actual_daily_deaths")) %>%
  pivot_wider(names_from = variable, values_from = mean_value) %>%
  rename(cases = actual_daily_cases,
         deaths = actual_daily_deaths) %>%
  mutate(deaths = ifelse(is.na(deaths), 0, deaths)) %>%
  bind_cols(hosps = NA_real_) %>%
  mutate(tmp = as.integer(date))

origin <- mdata2$tmp[which(mdata2$date == "2020-03-24")] - 21 #using 2020-03-03 - 0 would seem more natural, but that date is not in all data sets
mdata2$time <- mdata2$tmp - origin %>% as.integer()
mdata3 <- mdata2 %>% select(date, location, cases, hosps, deaths, time) %>% 
  filter(time >= 1) # time 1 is used as t0 in makepompmodel()
pdata_versioned <- mdata3

n_knots <- round(nrow(pdata_versioned) / 21)
knot_coefs <-  paste0("b", 1:n_knots)

max_obs_date <- max(pdata_versioned$date)

cdata <- mdata %>% filter(variable == "mobility_trend") %>% group_by(date) %>% 
  slice(1) %>% 
  mutate(time = as.integer(date) - as.integer(as.Date("2020-03-03"))) %>%
  filter(date <= max_obs_date) %>%
  filter(time >= 1)

covar <- covariate_table(
  t = pdata_versioned$time,
  seas = bspline.basis(
    x=t,
    nbasis=n_knots,
    degree=3
  ),
  rel_beta_change = as.matrix(cdata$mean_value),
  trend_sim = as.matrix(rep(10, times = nrow(cdata))),  # this is a placeholder only needed for simulation
  fit = 1,  # 1 = fitting; 0 = simulating
  times="t",
  order = "constant"
)

# --------------------------------------------------
# Specify parameters and initial state variables that were estimated
# --------------------------------------------------

est_these_pars = c("log_sigma_dw", "min_frac_dead", "max_frac_dead", "log_half_dead",
                   "log_theta_cases", "log_theta_deaths")
est_these_inivals = c("E1_0", "Ia1_0", "Isu1_0", "Isd1_0")

state_pops <- readRDS("../data/us_popsize.rds")
statedf <- state_pops %>% 
  # R0 at beginning of epidemic for each state
  dplyr::mutate(initR0 = dplyr::case_when(
    state_full %in% c("New York") ~ 10, 
    state_full %in% c("Illinois") ~ 8,
    state_full %in% c("Indiana") ~ 6, 
    state_full %in% c("Maryland") ~ 8,
    state_full %in% c("Massachusetts") ~ 6,
    state_full %in% c("New Jersey") ~ 6,
    state_full %in% c("Ohio") ~ 6,
    TRUE ~ 6 # default initial R0
  ))


# Set the parameter values and initial conditions
par_var_list <- setparsvars_warm(iniparvals = "fresh", # list or "fresh"
                                 est_these_pars = c(est_these_pars, knot_coefs), 
                                 est_these_inivals = est_these_inivals,
                                 population = statedf %>% 
                                   filter(state_full == locname) %>% pull(total_pop),
                                 n_knots = n_knots,
                                 # set R0 at beginning of epidemic
                                 rnaught = statedf %>% 
                                   filter(state_full == locname) %>% pull(initR0))



# Make the pomp model
pomp_res <- list(
  pomp_model = makepompmodel(
    par_var_list = par_var_list,
    pomp_data = pdata_versioned,
    pomp_covar = covar,
    n_knots = n_knots
  ),
  pomp_data = pdata_versioned,
  pomp_covar = covar,
  location = locname,
  par_var_list = par_var_list
)

# Run scenarios

keyvars <- c("MIF_ID", "LogLik", "LogLik_SE")
pars_path <- file.path("archive", fdt, paste0(locname, "-params-natural.rds"))
pn1 <-
  readRDS(pars_path)[keyvars,]
pn2 <- readRDS(pars_path) %>%
  filter(is_fitted == "yes")

assumed_fitted <-
  c(
    "min_frac_dead",
    "max_frac_dead",
    "log_half_dead",
    "theta_cases",
    "theta_deaths",
    "sigma_dw",
    "b1",
    "b2",
    "b3",
    "b4",
    "b5",
    "b6",
    "b7",
    "b8",
    "b9",
    "b10",
    "b11",
    "b12",
    "E1_0",
    "Ia1_0",
    "Isu1_0",
    "Isd1_0"
  )
stopifnot(isTRUE(setequal(rownames(pn2), assumed_fitted)))


pn <- bind_rows(pn1, pn2) %>% select(-is_fitted) %>% t() %>%
  as_tibble() %>%
  mutate(
    log_theta_cases = log(theta_cases),
    log_theta_deaths = log(theta_deaths),
    log_half_dead = log(log_half_dead),
    log_sigma_dw = log(sigma_dw),
    min_frac_dead = log(1 / min_frac_dead - 1),
    max_frac_dead = log(1 / max_frac_dead - 1),
    E1_0 = log(E1_0),
    Ia1_0 = log(Ia1_0),
    Isu1_0 = log(Isu1_0),
    Isd1_0 = log(Isd1_0)
  ) %>%
  select(-theta_cases,-theta_deaths,-sigma_dw)

mle_pars_path <- file.path("archive", fdt, paste0("parameter-estimates-", locname, ".rds"))
mle_params <-
  readRDS(mle_pars_path)
others <- setdiff(names(mle_params), colnames(pn))
pomp_res$all_partable <-
  bind_cols(pn, as_tibble(as.list(mle_params[others]))) %>%
  arrange(-LogLik)

pomp_res$scenarios <- runscenarios(
  pomp_res,
  par_var_list = pomp_res$par_var_list,
  forecast_horizon_days = 6 * 7
)
outdir <- file.path("weekly-forecast-simulations", fdt)
if(!dir.exists(outdir)){
  dir.create(outdir, recursive = TRUE)
}
filename <- file.path(outdir, paste0(locname, ".rds")) 
saveRDS(object = pomp_res, file = filename)
