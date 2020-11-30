library(evalcast)
library(magrittr)
library(ggplot2)
library(purrr)

# collect prediction cards
cpred <- list()
cpred$sqs <- load_local_covidhub("CEID-compart_mif_sq")
cpred$lis <- load_local_covidhub("CEID-compart_mif_li")
cpred$rns <- load_local_covidhub("CEID-compart_mif_rn")

fdts <- map(cpred[[1]], attr, "forecast_date") %>% 
  map_chr(as.character) %>% unique()

cpred$rwf <- get_covidhub_predictions("CEID-Walk",
                                      forecast_dates = fdts)
cpred$utf <- get_covidhub_predictions("UT-Mobility",
                                          forecast_dates = fdts)
cpred$cef <- get_covidhub_predictions("COVIDhub-ensemble",
                                          forecast_dates = fdts)

cpred2 <- map(cpred, filter_predictions, geo_type = "state")


# case forecast evaluations
pull_cases <- function(x, h = c(1, 2)){
  filter_predictions(x, 
                     response_signal = "confirmed_incidence_num", 
                     ahead = h)
}

cpred3 <- cpred2 %>% map(pull_cases)

case_evals <- cpred3 %>% map(evaluate_predictions, backfill_buffer = 0)

ci <- list()
ci$h1 <- evalcast:::intersect_locations(map(case_evals, 1))
ci$h2 <- evalcast:::intersect_locations(map(case_evals, 2))

add_ci <- function(p){
  # confidence interval for median calculated by `boxplot.stats`
  # idea from koshke at https://stackoverflow.com/a/8135865
  
  f <- function(x) {
    ans <- boxplot.stats(x)
    data.frame(ymin = ans$conf[1], ymax = ans$conf[2], y = ans$stat[3])
  }
  p + stat_summary(fun.data = f, geom = "crossbar", color = NA, 
                   fill = "skyblue", size = 5, alpha = 0.5, width = 0.75)
}

plot_measure(ci$h1, "ae") %>% add_ci()
plot_measure(ci$h2, "wis") %>% add_ci()


plot_measure(ci$h1, "wis") %>% add_ci()
plot_measure(ci$h2, "wis") %>% add_ci()
plot_width(ci$h1, levels = 0.9)
plot_width(ci$h2, levels = 0.9)

map(ci$h1, plot_calibration)


# death forecast evaluations

esqs2 <- evalcast::evaluate_predictions(sqs[2], backfill_buffer = 0)
elis2 <- evalcast::evaluate_predictions(lis[2], backfill_buffer = 0)
erns2 <- evalcast::evaluate_predictions(rns[2], backfill_buffer = 0)
erwf2 <- evaluate_predictions(srwf[2], backfill_buffer = 0)
eutf <- evaluate_predictions(utf[1], backfill_buffer = 0)
ecef2 <- evaluate_predictions(scef[2], backfill_buffer = 0)

ei2 <- evalcast:::intersect_locations(c(erwf2, esqs2, elis2, erns2, eutf, ecef2))
plot_measure(ei2, "ae") %>% add_ci()
plot_measure(ei2, "wis") %>% add_ci()
plot_width(ei2)
plot_calibration(ei2[[1]])
plot_calibration(ei2[[2]])
plot_calibration(ei2[[3]])
plot_calibration(ei2[[4]])
plot_calibration(ei2[[5]])
plot_calibration(ei2[[6]])
