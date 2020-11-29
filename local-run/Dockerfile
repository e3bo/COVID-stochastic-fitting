FROM docker.io/rocker/geospatial:4.0.2
MAINTAINER Eamon O'Dea <[last name without apostrophe]35@gmail.com>

RUN install2.r --skipinstalled --error \
doParallel \
dplyr \
ggplot2 \
here \
lubridate \
pomp \
pomp \
purrr \
readr \
tibble \
tidyr \
tidyverse \
vctrs

RUN installGithub.r e3bo/covidcast/R-packages/evalcast@f01331467bfc6c641a0db62a8c41fb7e5a3b5660
