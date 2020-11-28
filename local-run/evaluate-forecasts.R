library(evalcast)

sqs <- load_local_covidhub("CEID-compart_mif_sq")
rwf <- evalcast::get_covidhub_predictions("CEID-Walk", 
                                          forecast_dates = "2020-11-16")
srwf <- filter_predictions(rwf, geo_type = "state")

# case forecast evaluations
esqs <- evalcast::evaluate_predictions(sqs[1], backfill_buffer = 0)
erwf <- evaluate_predictions(srwf[1], backfill_buffer = 0)

ei <- evalcast:::intersect_locations(c(erwf, esqs))
plot_measure(ei, "ae")
plot_measure(ei, "wis")
plot_width(ei)
plot_calibration(ei[[1]])
plot_calibration(ei[[2]])

# death forecast evaluations

esqs2 <- evalcast::evaluate_predictions(sqs[2], backfill_buffer = 0)
erwf2 <- evaluate_predictions(srwf[2], backfill_buffer = 0)

ei2 <- evalcast:::intersect_locations(c(erwf2, esqs2))
plot_measure(ei2, "ae")
plot_measure(ei2, "wis")
plot_width(ei2)
plot_calibration(ei2[[1]])
plot_calibration(ei2[[2]])
