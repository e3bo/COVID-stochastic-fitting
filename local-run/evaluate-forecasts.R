library(dplyr)
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

fdtsrw <- c(paste0("2020-11-", c("16", "09", "01")),
            paste0("2020-10-", c("25", "18", "11", "04")),
            paste0("2020-09-", c("27", "20")))
cpred$rwf <- get_covidhub_predictions("CEID-Walk",
                                      forecast_dates = fdtsrw)
cpred$utf <- get_covidhub_predictions("UT-Mobility",
                                          forecast_dates = fdts)
cpred$cef <- get_covidhub_predictions("COVIDhub-ensemble",
                                          forecast_dates = fdts)

cpred2 <- map(cpred, filter_predictions, geo_type = "state")

## trajectory_plots

sqselc <- cpred$sqs %>% 
  filter_predictions(response_signal = "confirmed_incidence_num", 
                     ahead = c(1, 2))
rwselc <- cpred$rwf %>% 
  filter_predictions(response_signal = "confirmed_incidence_num", 
                     geo_type = "state", ahead = c(1, 2))

evalcast:::intersect_locations(c(sqselc, rwselc)) %>% 
  plot_trajectory(first_day = "2020-09-06")
ggsave("case-trajectories.png", width = 27, height = 18)

sqsel <- cpred$sqs %>% 
  filter_predictions(response_signal = "deaths_incidence_num", ahead = c(1, 2))
rwsel <- cpred$rwf %>% 
  filter_predictions(response_signal = "deaths_incidence_num", 
                     geo_type = "state", ahead = c(1, 2))

evalcast:::intersect_locations(c(sqsel, rwsel)) %>% 
  plot_trajectory(first_day = "2020-09-06")
ggsave("death-trajectories.png", width = 27, height = 18)

# case forecast evaluations
pull_cases <- function(x, h = c(1, 2)){
  filter_predictions(x, 
                     response_signal = "confirmed_incidence_num", 
                     ahead = h)
}

cpred3 <- cpred2[-5] %>% map(pull_cases)

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
                   fill = "skyblue", size = 5, alpha = 0.5, width = 0.75) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
}

plot_measure(ci$h1, "ae") %>% add_ci()
plot_measure(ci$h1, "wis") %>% add_ci()

plot_measure(ci$h2, "ae") %>% add_ci()
plot_measure(ci$h2, "wis") %>% add_ci()

my_dotplot <- function(sc, err_name = "wis"){
  nm <- attr(sc, "name_of_forecaster")
  rsp <- attr(sc, "response")$signal
  sc2 <- sc %>% mutate(location = covidcast::fips_to_abbr(location))
  ordered_levels <- sc2 %>% group_by(.data$location) %>% 
    summarize(avg_err = mean(.data[[err_name]])) %>% 
    arrange(.data$avg_err)
  sc2$location <- factor(sc2$location, levels = ordered_levels$location)

  sc2 %>% ggplot(aes(x = wis, y = location)) + geom_point(alpha = 0.5) + 
    scale_x_continuous(trans = "log1p") + 
    ggtitle(nm, subtitle = rsp) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
}

my_dotplot(ci$h1[[1]], "wis")
my_dotplot(ci$h1[[4]], "wis")
  
plot_width(ci$h1, levels = 0.9)
plot_width(ci$h2, levels = 0.9)

map(ci$h1, plot_calibration)

# death forecast evaluations

pull_deaths <- function(x, h = c(1, 2)){
  filter_predictions(x, 
                     response_signal = "deaths_incidence_num", 
                     ahead = h)
}

cpred4 <- cpred2 %>% map(pull_deaths)

death_evals <- cpred4 %>% map(evaluate_predictions, backfill_buffer = 1)

di <- list()
di$h1 <- evalcast:::intersect_locations(map(death_evals, 1))
di$h2 <- evalcast:::intersect_locations(map(death_evals, 2))

plot_measure(di$h1, "ae") %>% add_ci()
plot_measure(di$h1, "wis") %>% add_ci()

plot_measure(di$h2, "ae") %>% add_ci()
plot_measure(di$h2, "wis") %>% add_ci()

my_dotplot(di$h1[[1]], "wis")
my_dotplot(di$h1[[4]], "wis")

plot_width(di$h1, levels = 0.9)
plot_width(di$h2, levels = 0.9)

map(di$h1, plot_calibration)

## regression modeling

fit_mod <- function(sc){
  df <- dplyr::bind_rows(sc, .id = "forecaster")
  m <- lme4::lmer(log(wis + 1) ~ forecaster + (1|location) + (1|end), data = df)
  ci <- confint(m)
  list(model = m, conf = ci)  
}

(di_mods <- map(di, fit_mod))
(ci_mods <- map(ci, fit_mod))


