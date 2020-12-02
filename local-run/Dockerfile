FROM docker.io/rocker/geospatial:4.0.2
MAINTAINER Eamon O'Dea <[last name without apostrophe]35@gmail.com>

RUN install2.r --skipinstalled --error \
doParallel \
dplyr \
ggplot2 \
here \
lme4 \
lubridate \
pomp \
pomp \
purrr \
readr \
tibble \
tidyr \
tidyverse \
vctrs \
xml2

RUN installGithub.r e3bo/covidcast/R-packages/evalcast@93e56fd  
RUN apt-get update && \
    apt-get install -y --no-install-recommends --fix-missing \
        subversion
