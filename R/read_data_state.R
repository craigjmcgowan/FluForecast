# Pull in and save datasets from EpiData and other sources 
library(tidyverse)
library(cdcfluview)
library(MMWRweek)
library(gtrendsR)
library(lubridate)
library(FluSight)
library(epiforecast)

# Load functions
source("R/utils.R")
source('R/EpiDataAPI.R')

# Save current MMWR week in format Epidata wants
current_week <- as.numeric(paste0(MMWRweek(Sys.Date())[[1]], 
                                  str_pad(MMWRweek(Sys.Date())[[2]], width = 2,
                                          side = "left", pad = "0")))

pull_week <- case_when(
  MMWRweek(Sys.Date())[3] == 7 & substr(current_week, 5, 6) == "01" ~ 
    as.numeric(paste0(as.numeric(substr(current_week, 1, 4)) - 1, "52")),
  MMWRweek(Sys.Date())[3] == 7 ~
    current_week - 1,
  substr(current_week, 5, 6) == "01" ~ 
    as.numeric(paste0(as.numeric(substr(current_week, 1, 4)) - 1, "51")),
  substr(current_week, 5, 6) == "02" ~ 
    as.numeric(paste0(as.numeric(substr(current_week, 1, 4)) - 1, "52")),
  TRUE ~ current_week - 2
)

# Read in current state data
state_current <- Epidata$fluview(regions = list("pa", "de", "md",
                                             "va", "wv"),
                           epiweeks = list(Epidata$range(201040, pull_week)))$epidata %>%
  modify_depth(2, function(x) ifelse(is.null(x), NA, x)) %>%
  bind_rows() %>%
  mutate(year = as.integer(substr(epiweek, 1, 4)),
         week = as.integer(substr(epiweek, 5, 6)),
         season = ifelse(week >= 40,
                         paste0(year, "/", year + 1),
                         paste0(year - 1, "/", year)),
         location = case_when(
           region == "pa" ~ "Pennsylvania",
           region == "de" ~ "Delaware",
           region == "md" ~ "Maryland",
           region == "va" ~ "Virginia",
           region == "wv" ~ "West Virginia"
         )) %>%
  rename(ILI = ili) %>%
  # Make ILI always > 0
  mutate(ILI = case_when(near(ILI, 0) ~ 0.1,
                         TRUE ~ ILI)) %>%
  select(season, location, year, week, epiweek, issue, release_date, lag, ILI)

# Read in prior state data
epidata.cache.dir = "~/.epiforecast-cache"
if (!dir.exists(epidata.cache.dir)) {
  dir.create(epidata.cache.dir)
}

state_history = bind_rows(
  fetchEpidataHistoryDF(
    "fluview", "pa", 0:51,
    first.week.of.season = 40L,
    cache.file.prefix=file.path(epidata.cache.dir,paste0("fluview_pa"))
  ),
  fetchEpidataHistoryDF(
    "fluview", "de", 0:51,
    first.week.of.season = 40L,
    cache.file.prefix=file.path(epidata.cache.dir,paste0("fluview_de"))
  ),
  fetchEpidataHistoryDF(
    "fluview", "md", 0:51,
    first.week.of.season = 40L,
    cache.file.prefix=file.path(epidata.cache.dir,paste0("fluview_md"))
  ),
  fetchEpidataHistoryDF(
    "fluview", "va", 0:51,
    first.week.of.season = 40L,
    cache.file.prefix=file.path(epidata.cache.dir,paste0("fluview_va"))
  ),
  fetchEpidataHistoryDF(
    "fluview", "wv", 0:51,
    first.week.of.season = 40L,
    cache.file.prefix=file.path(epidata.cache.dir,paste0("fluview_wv"))
  )
) %>%
  filter(!is.na(ili)) %>%
  mutate(year = as.integer(substr(epiweek, 1, 4)),
         week = as.integer(substr(epiweek, 5, 6)),
         season = ifelse(week >= 40,
                         paste0(year, "/", year + 1),
                         paste0(year - 1, "/", year)),
         location = case_when(
           region == "pa" ~ "Pennsylvania",
           region == "de" ~ "Delaware",
           region == "md" ~ "Maryland",
           region == "va" ~ "Virginia",
           region == "wv" ~ "West Virginia"
         )) %>%
  rename(ILI = ili) %>%
  # Make ILI always > 0
  mutate(ILI = case_when(near(ILI, 0) ~ 0.1,
                         TRUE ~ ILI)) %>%
  select(season, location, year, week, epiweek, issue, release_date, lag, ILI)

# Create measure of backfill
state_backfill <- group_by(state_history, location, epiweek) %>%
  filter(year > 2017 | (year == 2017 & week >= 40)) %>%
  arrange(issue, .by_group = T) %>%
  summarize(backfill = last(ILI) - first(ILI))

save(state_current, state_history, state_backfill, 
     file = "Data/ili_state.Rdata")

# Virologic data ------
virologic_state <- who_nrevss(region = "state", years = c(2010:2018))

virologic_before_1516 <- virologic_state[[1]] %>%
  mutate_at(c("a_2009_h1n1", "a_h1", "a_h3", "a_subtyping_not_performed",
              "a_unable_to_subtype", "b", "h3n2v"),
            as.integer)

virologic_ph_lab <- virologic_state[[2]] %>%
  mutate_at(c("a_2009_h1n1", "a_h3", "a_subtyping_not_performed",
               "b", "bvic", "byam", "h3n2v"),
            as.integer)

virologic_state <- bind_rows(
  virologic_before_1516, virologic_ph_lab
) %>%
  # Create season indicator
  mutate(season = ifelse(week >= 40,
                         paste0(year, "/", year + 1),
                         paste0(year - 1, "/", year)),
         location = region) %>%
  filter(location %in% c("Pennsylvania", "Delaware", "Maryland", "Virginia",
                         "West Virginia")) %>%
  # Create mock H1 and H3 percents if all A samples had been tested
  rowwise() %>%
  mutate(pos_samples = sum(a_h1, a_2009_h1n1, a_h3, a_subtyping_not_performed, 
                           a_unable_to_subtype, b, bvic, byam, na.rm = TRUE),
         a_pos_samples = sum(a_h1, a_2009_h1n1, a_h3, a_subtyping_not_performed,
                             a_unable_to_subtype, na.rm = TRUE),
         b_pos_samples = sum(b, bvic, byam, na.rm = TRUE)) %>% #,
  # h1_per_samples = (sum(a_h1, a_2009_h1n1, na.rm = TRUE) / 
  #                     sum(a_h1, a_2009_h1n1, a_h3, na.rm = TRUE)) *
  #   sum(a_h1, a_2009_h1n1, a_h3, a_subtyping_not_performed, 
  #       a_unable_to_subtype, na.rm = T) / pos_samples,
  # h3_per_samples = a_h3 / sum(a_h1, a_2009_h1n1, a_h3, na.rm = TRUE) *
  #   sum(a_h1, a_2009_h1n1, a_h3, a_subtyping_not_performed, 
  #       a_unable_to_subtype, na.rm = T) / pos_samples,
  # b_per_samples = b_pos_samples / pos_samples) %>%
  ungroup() %>%
  mutate_at(vars(c("a_2009_h1n1", "a_h1", "a_h3", "a_subtyping_not_performed",
                   # "h1_per_samples", "h3_per_samples", "b_per_samples",
                   "a_unable_to_subtype", "b", "bvic", "byam")),
            function(x) ifelse(is.na(x), 0, x)) %>%
  # Create cumulative influenza percentage measures
  group_by(season, region) %>%
  arrange(wk_date, .by_group = T) %>%
  mutate(h1sum = cumsum(a_2009_h1n1) + cumsum(a_h1),
         h3sum = cumsum(a_h3),
         asum = cumsum(a_pos_samples),
         bsum = cumsum(b_pos_samples),
         # Subtype percentages of typed As
         h1per_of_a = ifelse(is.na(h1sum / (h1sum + h3sum)), 0.5, 
                             h1sum / (h1sum + h3sum)),
         h3per_of_a = ifelse(is.na(h3sum / (h1sum + h3sum)), 0.5, 
                             h3sum / (h1sum + h3sum)),
         # Cumulative percentages of each type, assuming A type %s are representative
         cum_bper = bsum / (asum + bsum),
         cum_h1per = h1per_of_a * asum / (asum + bsum),
         cum_h3per = h3per_of_a * asum / (asum + bsum)) %>%
  # If no samples reported, make each type 33%
  mutate(cum_h1per = ifelse(is.na(cum_h1per), 1/3, cum_h1per),
         cum_h3per = ifelse(is.na(cum_h3per), 1/3, cum_h3per),
         cum_bper = ifelse(is.na(cum_bper), 1/3, cum_bper)) %>% 
  ungroup() %>%
  select(location = region, season, year, week, a_h1, a_2009_h1n1, a_h3, 
         a_subtyping_not_performed, a_unable_to_subtype, b, bvic, byam, 
         h1per_of_a, h3per_of_a, cum_h1per, cum_h3per, cum_bper)

save(virologic_state,
     file = "Data/virologic_state.Rdata")