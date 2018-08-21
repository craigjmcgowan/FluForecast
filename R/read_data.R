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
                         paste0(year - 1, "/", year)),
         location = case_when(
           region == "nat" ~ "US National",
           region == "hhs1" ~ "HHS Region 1",
           region == "hhs2" ~ "HHS Region 2",
           region == "hhs3" ~ "HHS Region 3",
           region == "hhs4" ~ "HHS Region 4",
           region == "hhs5" ~ "HHS Region 5",
           region == "hhs6" ~ "HHS Region 6",
           region == "hhs7" ~ "HHS Region 7",
           region == "hhs8" ~ "HHS Region 8",
           region == "hhs9" ~ "HHS Region 9",
           region == "hhs10" ~ "HHS Region 10"
         )) %>%
  rename(ILI = wili) %>%
  select(-ili)

ili_current <- bind_rows(
  pull_curr_epidata(199740, 200139),
  pull_curr_epidata(200140, 200439),
  pull_curr_epidata(200440, 200730),
  pull_curr_epidata(200740, 201039),
  pull_curr_epidata(201040, 201339),
  pull_curr_epidata(201340, 201639),
  pull_curr_epidata(201640, this_week)) %>%
  mutate(year = as.integer(substr(epiweek, 1, 4)),
         week = as.integer(substr(epiweek, 5, 6)),
         season = ifelse(week >= 40,
                         paste0(year, "/", year + 1),
                         paste0(year - 1, "/", year)),
         location = case_when(
           region == "nat" ~ "US National",
           region == "hhs1" ~ "HHS Region 1",
             region == "hhs2" ~ "HHS Region 2",
             region == "hhs3" ~ "HHS Region 3",
             region == "hhs4" ~ "HHS Region 4",
             region == "hhs5" ~ "HHS Region 5",
             region == "hhs6" ~ "HHS Region 6",
             region == "hhs7" ~ "HHS Region 7",
             region == "hhs8" ~ "HHS Region 8",
             region == "hhs9" ~ "HHS Region 9",
             region == "hhs10" ~ "HHS Region 10"
         )) %>%
  rename(ILI = wili) %>%
  select(-ili)
  
save(ili_orig, ili_current, 
     file = "Data/ili.Rdata")

# Fetch virologic data from CDC -----
virologic_national <- who_nrevss(region = "national", years = c(1997:2017))
virologic_region <- who_nrevss(region = "hhs", years = c(1997:2017))

virologic_before_1516 <- bind_rows(
  virologic_national[[1]], virologic_region[[1]]
)

virologic_ph_lab <- bind_rows(
  virologic_national[[2]], virologic_region[[2]]
)

virologic_combined <- bind_rows(
  virologic_before_1516, virologic_ph_lab
) %>%
  # Create season indicator
  mutate(season = ifelse(week >= 40,
                         paste0(year, "/", year + 1),
                         paste0(year - 1, "/", year)),
         region = ifelse(region == "National", "US National", 
                         paste("HHS", region))) %>%
  rename(location = region)

save(virologic_combined,
     file = "Data/virologic.Rdata")

# Fetch Google Trends data -----
US_flu_0407 <- gtrends(keyword = "influenza",
                       geo = "US",
                       time = paste(MMWRweek2Date(2004, 41), MMWRweek2Date(2007, 40)))$interest_over_time

US_flu_0611 <- gtrends(keyword = "influenza",
                       geo = "US",
                       time = paste(MMWRweek2Date(2006, 41), MMWRweek2Date(2011, 40)))$interest_over_time

US_flu_1015 <- gtrends(keyword = "influenza",
                      geo = "US",
                      time = paste(MMWRweek2Date(2010, 41), MMWRweek2Date(2015, 40)))$interest_over_time

US_flu_1419 <- gtrends(keyword = "influenza",
                       geo = "US")$interest_over_time
gdata1 <- US_flu_0407
gdata2 <- US_flu_0611
US_flu_ratio <- function(gdata1, gdata2) {
  inner_join(gdata1 %>%
               select(date, old_hits = hits) %>%
               filter(MMWRweek(date)[[2]] > 40 | MMWRweek(date)[[2]] < 20),
             gdata2 %>%
               select(date, new_hits = hits) %>%
               filter(MMWRweek(date)[[2]] > 40 | MMWRweek(date)[[2]] < 20),
             by = "date") %>%
  mutate(ratio = old_hits / new_hits) %>%
  pull(ratio) %>%
  mean()
}

gratio_0607 <- US_flu_ratio(US_flu_0407, US_flu_0611)
gratio_1011 <- US_flu_ratio(US_flu_0611, US_flu_1015)
gratio_1415 <- US_flu_ratio(US_flu_1015, US_flu_1419)

gtrend_US_flu_merge <- filter(US_flu_0407, !date %in% US_flu_0611$date) %>%
  mutate(hits = hits / gratio_0607) %>%
  bind_rows(US_flu_0611) %>%
  filter(!date %in% US_flu_1015$date) %>%
  mutate(hits = hits / gratio_1011) %>%
  bind_rows(US_flu_1015) %>%
  filter(!date %in% US_flu_1419$date) %>%
  mutate(hits = hits / gratio_1415) %>%
  bind_rows(US_flu_1419) %>%
  select(date, hits) %>%
  mutate(year = MMWRweek(date)[[1]],
         week = MMWRweek(date)[[2]],
         season = ifelse(week >= 40,
                         paste0(year, "/", year + 1),
                         paste0(year - 1, "/", year)))

save(gtrend_US_flu_merge, file = "Data/Gtrends.Rdata")
