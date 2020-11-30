library(evalcast)
library(magrittr)
library(ggplot2)

sqs <- load_local_covidhub("CEID-compart_mif_sq")
lis <- load_local_covidhub("CEID-compart_mif_li")
rns <- load_local_covidhub("CEID-compart_mif_rn")

rwf <- evalcast::get_covidhub_predictions("CEID-Walk", 
                                          forecast_dates = "2020-11-16")
utf <- evalcast::get_covidhub_predictions("UT-Mobility",
                                          forecast_dates = "2020-11-16")
cef <- evalcast::get_covidhub_predictions("COVIDhub-ensemble",
                                          forecast_dates = "2020-11-16")
scef <- filter_predictions(cef, geo_type = "state")
srwf <- filter_predictions(rwf, geo_type = "state")

# case forecast evaluations
esqs <- evalcast::evaluate_predictions(sqs[1], backfill_buffer = 0)
elis <- evalcast::evaluate_predictions(lis[1], backfill_buffer = 0)
erns <- evalcast::evaluate_predictions(rns[1], backfill_buffer = 0)
ecef <- evaluate_predictions(scef[1], backfill_buffer = 0)
erwf <- evaluate_predictions(srwf[1], backfill_buffer = 0)


ei <- evalcast:::intersect_locations(c(erwf, esqs, elis, erns, ecef))

add_ci <- function(p){
  # confidence interval for median calculated by `boxplot.stats`, idea from koshke at https://stackoverflow.com/a/8135865
  
  f <- function(x) {
    ans <- boxplot.stats(x)
    data.frame(ymin = ans$conf[1], ymax = ans$conf[2], y = ans$stat[3])
  }
  p + stat_summary(fun.data = f, geom = "crossbar", color = NA, fill = "skyblue", size = 5, alpha = 0.5, width = 0.75)
}
plot_measure(ei, "ae") %>% add_ci()
plot_measure(ei, "wis") %>% add_ci()
plot_width(ei)
plot_calibration(ei[[1]])
plot_calibration(ei[[2]])
plot_calibration(ei[[3]])
plot_calibration(ei[[4]])
plot_calibration(ei[[5]])

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
