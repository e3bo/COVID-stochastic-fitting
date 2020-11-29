#!/usr/bin/env Rscript

library(tidyverse)

quant <- function(x, p){
  quantile(x, prob = p, names = FALSE, type = 8, na.rm = TRUE)
}

vardf <- function(var, samp){
  cname <- switch(var,
                  "inc hosp" = "hosps",
                  "inc death" = "deaths",
                  "inc case" = "cases",
                  "cum death" = "cum_deaths")
  if (var != "inc case") {
    prob <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)
  } else {
    prob <- c(0.025, 0.100, 0.250, 0.500, 0.750, 0.900, 0.975)
  }
  t1 <- tibble(
    target = var,
    type  = "quantile",
    quantile = prob,
    value = quant(samp[[cname]], prob)
  )
  t2 <-
    tibble(
      target = var,
      type = "point",
      quantile = NA,
      value = quant(samp[[cname]], 0.5)
      
    )
  bind_rows(t1, t2)
}

samp_to_df <-
  function(sampdf,
           vars = c("inc hosp", "inc death", "cum death", "inc case")) {
    purrr::map_dfr(vars, vardf, samp = sampdf)
  }

# Take simulation trajectories and output a data frame in the format described
# here: https://github.com/reichlab/covid19-forecast-hub/blob/6a7e5624ef540a55902770b7c17609d19e1f593a/data-processed/README.md
paths_to_forecast <- function(out, loc = "13", wks_ahead = 1:4, fdt) {
  if(any(wks_ahead > 20)){
    stop("Max weeks ahead accepted is 20", .call = FALSE)
  }
  
  out2 <- 
    out %>% group_by(Rep) %>% arrange(Date)
  
  take_back_step <- lubridate::wday(fdt, label = TRUE) %in% c("Sun", "Mon")
  week0 <- lubridate::epiweek(fdt) - take_back_step
  forecast_epiweeks <- week0 + wks_ahead
  
  weekly <- out2 %>%
    as_tibble() %>%
    mutate(epiweek = lubridate::epiweek(Date)) %>%
    filter(epiweek %in% forecast_epiweeks) %>%
    group_by(epiweek, Rep) %>%
    summarize(cases = sum(cases),
              deaths = sum(deaths),
              n = length(unique(Date)),
              Date = max(Date),
              .groups = "drop") %>%
    filter(n == 7) %>%
    rename("target_end_date" = Date) %>%
    nest(data = c(Rep, cases, deaths)) %>%
    mutate(pred_df = purrr::map(data, samp_to_df, 
                                vars = c("inc death", "inc case"))) %>%
    select(-data) %>%
    unnest(pred_df) %>%
    mutate(target = paste(epiweek - week0, 
                          "wk ahead", target)) %>%
    add_column(location = loc) %>%
    mutate(quantile = round(quantile, digits = 3),
           value = round(value, digits = 3)) %>%
    add_column(forecast_date = fdt) %>%
    select(forecast_date, target, target_end_date, location, type, quantile, 
           value)
  
  weekly %>% 
    filter(!is.na(value)) %>%
    filter((nchar(location)) <= 2 | str_detect(target, "inc case$")) ## only case forecasts accepted for counties
}

file2df <- function(fn, fdt, stp){
  res <- readRDS(fn)
  sim <- res$scenarios$sims
  sim_start <- sim %>% filter(Period == "Future") %>% pull(Date) %>% min()
  stopifnot(as.Date(fdt) >= sim_start)
  simout <- sim %>% 
    select(SimType, Date, Rep = rep_id, cases = C_new, deaths = D_new)
  sq <- filter(simout, SimType == stp) %>% select(-SimType)
  fcst <- paths_to_forecast(sq, fdt = fdt)
  fcst
}

write_scenario_fcst <- function(fdt, stp){
  data_dir <- file.path("weekly-forecast-simulations", fdt)
  files <- dir(data_dir, full.names = TRUE)
  fcst <- map(files, file2df, fdt = fdt, stp = stp) %>% bind_rows()
  suffix <- switch(stp,
                   status_quo = "sq",
                   linear_increase_sd = "li",
                   return_normal = "rn")
  fname <- paste0(fdt, "-CEID-compart_mif_", suffix, ".csv")
  dirname <- paste0("CEID-compart_mif_", suffix)
  if(!dir.exists(dirname)){
    dir.create(dirname)
  }
  path <- file.path(dirname, fname)
  write_csv(fcst, path = path)
}

fdt <- "2020-11-16"
write_scenario_fcst(fdt, stp = "status_quo")
write_scenario_fcst(fdt, stp = "linear_increase_sd")
write_scenario_fcst(fdt, stp = "return_normal")
