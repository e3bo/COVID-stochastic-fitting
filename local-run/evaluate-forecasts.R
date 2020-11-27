load_local_covidhub <- function (covid_hub_forecaster_name, forecast_dates = NULL, ...) 
{
  require(tidyverse)
  # modification of evalcast::get_covidhub_predictions() to get predictions from locally stored csv
  pcards <- list()
  if (is.null(forecast_dates)){
    out <- dir(covid_hub_forecaster_name) %>% 
      stringr::str_match_all(sprintf("(20\\d{2}-\\d{2}-\\d{2})-%s.csv", 
                                   covid_hub_forecaster_name))
    forecast_dates <- lubridate::as_date(out[[1]][, 2])
  }
  forecast_dates <- as.character(forecast_dates)
  for (forecast_date in forecast_dates) {
    filename <- sprintf("%s/%s-%s.csv", covid_hub_forecaster_name, 
                        forecast_date, covid_hub_forecaster_name)
    message("Downloading ", filename)
    pred <-
      readr::read_csv(
        filename,
        col_types = readr::cols(
          location = readr::col_character(),
          forecast_date = readr::col_date(format = ""),
          quantile = readr::col_double(),
          value = readr::col_double(),
          target = readr::col_character(),
          target_end_date = readr::col_date(format = ""),
          type = readr::col_character()
        )
      )
    pcards[[forecast_date]] <- pred %>% rename(probs = .data$quantile, 
                                               quantiles = .data$value) %>% filter(str_detect(.data$target, 
                                                                                              "wk ahead inc")) %>% filter(.data$type == "quantile") %>% 
      separate(.data$target, into = c("ahead", NA, NA, 
                                      NA, "response"), remove = TRUE) %>% select(-.data$forecast_date, 
                                                                                 -.data$type, -.data$target_end_date) %>% mutate(geo_type = case_when(nchar(.data$location) == 
                                                                                                                                                        2 ~ "state", nchar(.data$location) == 5 ~ "county")) %>% 
      group_by(.data$ahead, .data$response, .data$location, 
               .data$geo_type) %>% group_modify(~tibble(forecast_distribution = list(.))) %>% 
      group_by(.data$ahead, .data$response, .data$geo_type) %>% 
      group_map(~{
        if (.y$response == "death") {
          signals <- tibble(data_source = "jhu-csse", 
                            signal = "deaths_incidence_num")
        }
        else if (.y$response == "case") {
          signals <- tibble(data_source = "jhu-csse", 
                            signal = "confirmed_incidence_num")
        }
        attributes(.x) <- c(attributes(.x), list(forecaster = NA, 
                                                 name_of_forecaster = covid_hub_forecaster_name, 
                                                 signals = signals, forecast_date = lubridate::ymd(forecast_date), 
                                                 incidence_period = "epiweek", ahead = as.numeric(.y$ahead), 
                                                 geo_type = .y$geo_type, geo_values = NA, from_covidhub = FALSE))
        class(.x) <- c("prediction_card", class(.x))
        return(.x)
      })
  }
  purrr::flatten(pcards) %>% evalcast::filter_predictions(...)
}

library(evalcast)

sqs <- load_local_covidhub("CEID-compart_mif_sq")
esqs <- evalcast::evaluate_predictions(sqs[2], backfill_buffer = 0)

rwf <- evalcast::get_covidhub_predictions("CEID-Walk", 
                                          forecast_dates = "2020-11-16")
srwf <- filter_predictions(rwf, geo_type = "state")
erwf <- evaluate_predictions(srwf[2], backfill_buffer = 0)

ei <- evalcast:::intersect_locations(c(erwf, esqs))


