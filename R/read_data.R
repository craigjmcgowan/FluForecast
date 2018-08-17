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

pull_curr_epidata <- function(start, end) {
  Epidata$fluview(list('nat', 'hhs1', 'hhs2', 'hhs3', 'hhs4', 'hhs5',
                       'hhs6', 'hhs7', 'hhs8', 'hhs9', 'hhs10'),
                  list(Epidata$range(start, end)))$epidata %>%
    modify_depth(2, function(x) ifelse(is.null(x), NA, x)) %>%
    bind_rows()
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
  bind_rows() %>%
  mutate(year = as.integer(substr(epiweek, 1, 4)),
         week = as.integer(substr(epiweek, 5, 6)),
         season = ifelse(week >= 40,
                         paste0(year, "/", year + 1),
                         paste0(year - 1, "/", year)))

ili_current <- pull_curr_epidata(201040, 201339) %>%
  bind_rows(pull_curr_epidata(201340, 201639)) %>%
  bind_rows(pull_curr_epidata(201640, this_week)) %>%
  mutate(year = as.integer(substr(epiweek, 1, 4)),
         week = as.integer(substr(epiweek, 5, 6)),
         season = ifelse(week >= 40,
                         paste0(year, "/", year + 1),
                         paste0(year - 1, "/", year)))
  
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

virologic_combined <- bind_rows(
  virologic_before_1516, virologic_ph_lab
) %>%
  # Create season indicator
  mutate(season = ifelse(week >= 40,
                         paste0(year, "/", year + 1),
                         paste0(year - 1, "/", year)),
         region = case_when(
           region == "National" ~ "nat",
           region == "Region 1" ~ "hhs1",
           region == "Region 2" ~ "hhs2",
           region == "Region 3" ~ "hhs3",
           region == "Region 4" ~ "hhs4",
           region == "Region 5" ~ "hhs5",
           region == "Region 6" ~ "hhs6",
           region == "Region 7" ~ "hhs7",
           region == "Region 8" ~ "hhs8",
           region == "Region 9" ~ "hhs9",
           region == "Region 10" ~ "hhs10"
         ))

save(virologic_combined,
     file = "Data/virologic.Rdata")

# Fetch Google Trends data -----
US_flu_1015 <- gtrends(keyword = "influenza",
                      geo = "US",
                      time = paste(MMWRweek2Date(2010, 41), MMWRweek2Date(2015, 40)))$interest_over_time

US_flu_1419 <- gtrends(keyword = "influenza",
                       geo = "US")$interest_over_time

US_flu_ratio <- inner_join(US_flu_1015 %>%
                             select(date, old_hits = hits) %>%
                             filter(date %within% interval(MMWRweek2Date(2013, 42), 
                                                           MMWRweek2Date(2014, 18)) |
                                      date %within% interval(MMWRweek2Date(2014, 42), 
                                                             MMWRweek2Date(2015, 18))),
                           US_flu_1419 %>%
                             select(date, new_hits = hits) %>%
                             filter(date %within% interval(MMWRweek2Date(2013, 42), 
                                                           MMWRweek2Date(2014, 18)) |
                                      date %within% interval(MMWRweek2Date(2014, 42), 
                                                             MMWRweek2Date(2015, 18))),
                           by = "date") %>%
  mutate(ratio = old_hits / new_hits) %>%
  pull(ratio) %>%
  mean()

US_flu_merge <- bind_rows(filter(US_flu_1015, !date %in% US_flu_1419$date),
                          mutate(US_flu_1419, hits = hits * US_flu_ratio)) %>%
  select(date, hits) %>%
  mutate(year = MMWRweek(date)[[1]],
         week = MMWRweek(date)[[2]],
         season = ifelse(week >= 40,
                         paste0(year, "/", year + 1),
                         paste0(year - 1, "/", year)))

save(US_flu_merge, file = "Data/Gtrends.Rdata")
