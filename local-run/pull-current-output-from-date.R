#!/usr/bin/env Rscript

library(magrittr)

fdt <- Sys.getenv("fdt")
fdtm <- as.POSIXct(paste(fdt, "0:03:00"), tz = "EST")

xmlfile <- tempfile()
url <- paste0("https://github.com/CEIDatUGA/COVID-stochastic-fitting.git",
              "/trunk/output/current")
cmd <- paste("svn log ", url, " --xml >", xmlfile)
system(cmd)
log <- xml2::read_xml(xmlfile)
datetimes <- 
  xml2::xml_find_all(log, "//logentry/date") %>% 
  xml2::xml_text() %>%
  lubridate::as_datetime()
is_later <- datetimes > fdtm
contemp <- match(FALSE, is_later)
xpath <- paste0("//logentry[", contemp, "]")
revision <- 
  xml2::xml_find_all(log, xpath) %>% 
  xml2::xml_attr("revision")
url2 <- paste0(url, "@r", revision)
system("mkdir archive")
cmd2 <- paste0("svn export ", url2, " ./archive/", fdt)
system(cmd2)
