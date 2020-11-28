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

source("./code/forward-simulations/simulate_trajectories.R")
source("./code/forward-simulations/runscenarios.R") #runs forward simulations 
source("./code/model-setup/makepompmodel.R") #generates the pomp model

# Load the pomp information
pomp_listr <- readRDS("./header/pomp_list.rds")
myargument <- 38
this_pomp <- pomp_listr[[myargument]]
fdt <- "2020-11-16"

locname <- this_pomp$location

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

origin <- mdata2$tmp[which(mdata2$date == "2020-03-04")] - 1
mdata2$time <- mdata2$tmp - origin %>% as.integer()
mdata3 <- mdata2 %>% select(date, location, cases, hosps, deaths, time)

# uncomment to test that data has not changed
#pdata_windowed <- this_pomp$pomp_data %>% filter(date < fdt)
#B <- mdata3 %>% filter(date > "2020-03-06")
#A <- pdata_windowed %>% filter(date < "2020-11-11")
#all.equal(A %>% select(-hosps), B %>% select(-hosps)) # hops are not used in fitting so exclude

pdata_versioned <- mdata3

covar_dates <- as.Date(this_pomp$pomp_covar@times, origin = "2020-03-03")
covar_windowed <- this_pomp$pomp_covar
covar_windowed@table <- covar_windowed@table[, covar_dates < fdt]
covar_windowed@times <- covar_windowed@times[covar_dates < fdt]

this_pomp$pomp_data <- pdata_versioned
this_pomp$pomp_covar <- covar_windowed

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
  readRDS("./archive/2020-11-16/Georgia-params-natural.rds")[keyvars,]
pn2 <- readRDS("./archive/2020-11-16/Georgia-params-natural.rds") %>%
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
  readRDS("./archive/2020-11-16/parameter-estimates-Georgia.rds")
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
