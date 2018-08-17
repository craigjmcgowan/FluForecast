# Pull in and save datasets from EpiData and other sources 
library(tidyverse)
library(cdcfluview)
library(MMWRweek)
library(gtrendsR)
library(lubridate)

# Helper function to turn null into NA
null_to_na <- function(x) {
  if (is.null(x)) x <- NA
}

# Save current MMWR week in format Epidata wants
this_week <- as.numeric(paste0(MMWRweek(Sys.Date())[[1]], MMWRweek(Sys.Date())[[2]] ))



# Fetch wILI data from EpiData API -------
source('R/EpiDataAPI.R')
ili_orig <- Epidata$fluview(list('nat', 'hhs1', 'hhs2', 'hhs3', 'hhs4', 'hhs5',
                                 'hhs6', 'hhs7', 'hhs8', 'hhs9', 'hhs10'),
                            list(Epidata$range(201401, this_week)),
                            lag = 0)$epidata %>%
  modify_depth(2, function(x) ifelse(is.null(x), NA, x)) %>%
  bind_rows()

ili_current <- Epidata$fluview(list('nat', 'hhs1', 'hhs2', 'hhs3', 'hhs4', 'hhs5',
                                    'hhs6', 'hhs7', 'hhs8', 'hhs9', 'hhs10'),
                               list(Epidata$range(201401, this_week)))$epidata %>%
  modify_depth(2, function(x) ifelse(is.null(x), NA, x)) %>%
  bind_rows()

save(ili_orig, ili_current, 
     file = "Data/ili.Rdata")

# Fetch virologic data from CDC -----
virologic_national <- who_nrevss(region = "national", years = c(2010:2017))
virologic_region <- who_nrevss(region = "hhs", years = c(2010:2017))

virologic_before_1516 <- bind_rows(
  virologic_national[[1]], virologic_region[[1]]
)

virologic_ph_lab <- bind_rows(
  virologic_national[[2]], virologic_region[[1]]
)

save(virologic_before_1516, virologic_ph_lab,
     file = "Data/virologic.Rdata")

# Fetch Google Trends data -----
US_flu_1015 <- gtrends(keyword = "influenza",
                      geo = "US",
                      time = paste(MMWRweek2Date(2010, 40), MMWRweek2Date(2015, 40)))$interest_over_time

US_flu_1419 <- gtrends(keyword = "influenza",
                       geo = "US")$interest_over_time

US_flu_ratio <- inner_join(US_flu_1015 %>%
                             select(date, old_hits = hits) %>%
                             filter(date %within% interval(MMWRweek2Date(2013, 42), 
                                                           MMWRweek2Date(2014, 18)) |
                                      date %within% interval(MMWRweek2Date(2014, 42), 
                                                             MMWRweek2Date(2015, 18))),
                           US_flu_1418 %>%
                             select(date, new_hits = hits) %>%
                             filter(date %within% interval(MMWRweek2Date(2013, 42), 
                                                           MMWRweek2Date(2014, 18)) |
                                      date %within% interval(MMWRweek2Date(2014, 42), 
                                                             MMWRweek2Date(2015, 18))),
                           by = "date") %>%
  mutate(ratio = old_hits / new_hits) %>%
  pull(ratio) %>%
  mean()

US_flu_merge <- bind_rows(filter(US_flu_1014, !date %in% US_flu_1418$date),
                          mutate(US_flu_1418, hits = hits * US_flu_ratio)) %>%
  select(date, hits)

save(US_flu_merge, file = "Data/Gtrends.Rdata")
