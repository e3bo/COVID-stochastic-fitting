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


# --------------------------------------------------
# Source all needed functions/scripts
# --------------------------------------------------

source("./code/forward-simulations/simulate_trajectories.R")
source("./code/forward-simulations/runscenarios.R") #runs forward simulations 
source("./code/model-setup/makepompmodel.R") #generates the pomp model

# Load the pomp information
pomp_listr <- readRDS("./header/pomp_list.rds")
myargument <- 38 ## Georgia
this_pomp <- pomp_listr[[myargument]]
n_knots <- round(nrow(this_pomp$pomp_data) / 21)


# --------------------------------------------------
# Create a time-stamp variable
# Will be applied to saved results
# defined outside loop so it's the same for each state
# --------------------------------------------------
timestamp <- readRDS("./header/timestamp.rds")

# Make the pomp model
pomp_model <- makepompmodel(
  par_var_list = this_pomp$par_var_list,
  pomp_data = this_pomp$pomp_data,
  pomp_covar = this_pomp$pomp_covar,
  n_knots = n_knots
)
this_pomp$pomp_model <- pomp_model

pomp_res <- this_pomp #current state
rm(this_pomp) #remove the old object

# Run scenarios

keyvars <- c("MIF_ID", "LogLik", "LogLik_SE")
pn1 <-
  readRDS("./output/current/Georgia-params-natural.rds")[keyvars,]
pn2 <- readRDS("./output/current/Georgia-params-natural.rds") %>%
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


mle_params <-
  readRDS("./output/current/parameter-estimates-Georgia.rds")
others <- setdiff(names(mle_params), colnames(pn))
pomp_res$all_partable <-
  bind_cols(pn, as_tibble(as.list(mle_params[others]))) %>%
  arrange(-LogLik)

pomp_res$scenarios <- runscenarios(
  pomp_res,
  par_var_list = pomp_res$par_var_list,
  forecast_horizon_days = 6 * 7
)
filename = paste0('./output/', pomp_res$filename_label, '_weekly_results.rds')
saveRDS(object = pomp_res, file = filename)
