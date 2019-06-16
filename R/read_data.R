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

state_matchup <- tibble(abb = tolower(state.abb),
                        name = state.name) %>%
  # Add matchups for regions and national
  bind_rows(
    tibble(abb = c("nat", "hhs1", "hhs2", "hhs3", "hhs4", "hhs5", "hhs6", 
                   "hhs7", "hhs8", "hhs9", "hhs10"),
           name = c("US National", "HHS Region 1", "HHS Region 2", "HHS Region 3",
                    "HHS Region 4", "HHS Region 5", "HHS Region 6",  "HHS Region 7", 
                    "HHS Region 8", "HHS Region 9", "HHS Region 10"))
  )


# ILI data ------
current_epidata <- function(start_wk, end_wk) {
  Epidata$fluview(
    regions = list("al", "ak", "az", "ar", "ca", "co", "ct", "de", "ga", "hi", 
                   "id", "il", "in", "ia", "ks", "ky", "la", "me", "md", "ma", 
                   "mi", "mn", "ms", "mo", "mt", "ne", "nv", "nh", "nj", "nm", 
                   "ny", "nc", "nd", "oh", "ok", "or", "pa", "ri", "sc", "sd", 
                   "tn", "tx", "ut", "vt", "va", "wa", "wv", "wi", "wy", "nat",
                   "hhs1", "hhs2", "hh3", "hhs4", "hhs5", "hhs6", "hhs7", "hhs8",
                   "hhs9", "hhs10"),
    epiweeks = list(Epidata$range(start_wk, end_wk))
  )$epidata %>%
    modify_depth(2, function(x) ifelse(is.null(x), NA, x)) %>%
    bind_rows()
}

ili_current <- bind_rows(
  current_epidata(200040, 200139),
  current_epidata(200140, 200239),
  current_epidata(200240, 200339),
  current_epidata(200340, 200439),
  current_epidata(200440, 200539),
  current_epidata(200540, 200639),
  current_epidata(200640, 200739),
  current_epidata(200740, 200839),
  current_epidata(200840, 200939),
  current_epidata(200940, 201039),
  current_epidata(201040, 201139),
  current_epidata(201140, 201239),
  current_epidata(201240, 201339),
  current_epidata(201340, 201439),
  current_epidata(201440, 201539),
  current_epidata(201540, 201639),
  current_epidata(201640, 201739),
  current_epidata(201740, 201839),
  current_epidata(201840, pull_week)
) %>%
  # Add in location name
  left_join(rename(state_matchup, location = name), by = c("region" = "abb")) %>%
  # Fix formatting of some variables
  mutate(year = as.integer(substr(epiweek, 1, 4)),
         week = as.integer(substr(epiweek, 5, 6)),
         season = ifelse(week >= 40,
                         paste0(year, "/", year + 1),
                         paste0(year - 1, "/", year)),
         ILI = case_when(region %in% state.abb ~ ili,
                         TRUE ~ wili),
         # Make ILI always > 0
         ILI = case_when(near(ILI, 0) ~ 0.1,
                         TRUE ~ ILI)) %>%
  select(season, location, year, week, epiweek, issue, release_date, lag, ILI)

# Read in prior state data
epidata.cache.dir = "~/.epiforecast-cache"
if (!dir.exists(epidata.cache.dir)) {
  dir.create(epidata.cache.dir)
}

epidata_history <- function(abb) {
  fetchEpidataHistoryDF(
    "fluview", abb, 0:51,
    first.week.of.season = 40L,
    cache.file.prefix=file.path(epidata.cache.dir,paste0(paste0("fluview_", abb)))
  )
}

history <- lapply(state_matchup$abb[state_matchup$abb != "fl"], epidata_history) %>%
  bind_rows() %>%
  distinct() %>%
  filter(!is.na(ili)) %>%
  # Add in state name
  left_join(rename(state_matchup, location = name), by = c("region" = "abb")) %>%
  # Fix formatting of some variables
  mutate(year = as.integer(substr(epiweek, 1, 4)),
         week = as.integer(substr(epiweek, 5, 6)),
         season = ifelse(week >= 40,
                         paste0(year, "/", year + 1),
                         paste0(year - 1, "/", year)),
         ILI = case_when(region %in% state.abb ~ ili,
                         TRUE ~ wili),
         # Make ILI always > 0
         ILI = case_when(near(ILI, 0) ~ 0.1,
                         TRUE ~ ILI)) %>%
  filter(year >= 2011 | season == "2010/2011") %>%
  select(season, location, year, week, epiweek, issue, release_date, lag, ILI)

# Create list of initial published values for coordination with modeling code
ili_init_pub_list <- list()

for(i in c(201040:201052, 201101:201152, 201201:201252, 201301:201338,
           201340:201352, 201401:201453, 201501:201552, 201601:201652,
           201701:201752, 201801:201852, 201901:pull_week)) {
  
  if(i %in% history$issue) {
    ili_init_pub_list[[paste(i)]] <- filter(history, issue == i)
  } else {
    ili_init_pub_list[[paste(i)]] <- history %>%
      filter(epiweek <= i, epiweek > i-100, issue == 201740)
  }
  
}

ili_orig <- filter(history, lag == 0)

# Create measure of backfill
ili_backfill <- left_join(
  select(ili_orig, location, season, year, week, orig_ILI = ILI),
  select(ili_current, location, season, year, week, final_ILI = ILI),
  by = c("location", "season", "year", "week")
) %>%
  mutate(backfill = final_ILI - orig_ILI)

ili_backfill_avg <- ili_backfill %>%
  group_by(location, week) %>%
  summarize(avg_backfill = mean(backfill))

saveRDS(ili_orig, "Data/ili_orig.RDS")
saveRDS(ili_current, "Data/ili_current.RDS")
saveRDS(ili_init_pub_list, "Data/ili_init_pub_list.RDS")
saveRDS(ili_backfill, "Data/ili_backfill.RDS")
saveRDS(ili_backfill_avg, "Data/ili_backfill_avg.RDS")

# Fetch virologic data from CDC -----
virologic_national <- who_nrevss(region = "national", years = c(1997:2018))
virologic_region <- who_nrevss(region = "hhs", years = c(1997:2018))
virologic_state <- who_nrevss(region = "state", years = c(2010:2018))

virologic_before_1516 <- bind_rows(
  virologic_national[[1]],
  virologic_region[[1]], 
  virologic_state[[1]] %>%
    mutate_at(c("a_2009_h1n1", "a_h1", "a_h3", "a_subtyping_not_performed",
                "a_unable_to_subtype", "b", "h3n2v", "total_specimens",
                "percent_positive"),
              as.integer)
)

virologic_ph_lab <- bind_rows(
  virologic_national[[2]], 
  virologic_region[[2]], 
  virologic_state[[2]] %>%
    mutate_at(c("a_2009_h1n1", "a_h3", "a_subtyping_not_performed",
                 "b", "bvic", "byam", "h3n2v", "total_specimens"),
              as.integer)
)

virologic_combined <- bind_rows(
  virologic_before_1516, virologic_ph_lab
) %>%
  # Create season indicator
  mutate(season = ifelse(week >= 40,
                         paste0(year, "/", year + 1),
                         paste0(year - 1, "/", year)),
         region = case_when(
           region == "National" ~ "US National", 
           str_detect(region, "Region") ~ paste("HHS", region),
           TRUE ~ region
         )) %>%
  # Create mock H1 and H3 percents if all A samples had been tested
  rowwise() %>%
  mutate(pos_samples = sum(a_h1, a_2009_h1n1, a_h3, a_subtyping_not_performed, 
           a_unable_to_subtype, b, bvic, byam, na.rm = TRUE),
         a_pos_samples = sum(a_h1, a_2009_h1n1, a_h3, a_subtyping_not_performed,
                             a_unable_to_subtype, na.rm = TRUE),
         b_pos_samples = sum(b, bvic, byam, na.rm = TRUE)) %>% 
  ungroup() %>%
  mutate_at(vars(c("a_2009_h1n1", "a_h1", "a_h3", "a_subtyping_not_performed",
                   "a_unable_to_subtype", "b", "bvic", "byam")),
            function(x) ifelse(is.na(x), 0, x)) %>%
  # Create cumulative influenza percentage measures
  group_by(season, region) %>%
  arrange(wk_date, .by_group = TRUE) %>%
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

saveRDS(virologic_combined, file = "Data/virologic.Rds")

# Fetch Google Trends data -----
gtrend_US_flu_merge <- fetch_gtrend("US") %>%
  mutate(location = "US National")

gtrend_state_list <- tibble(state_abb = state.abb) %>%
  mutate(data = map(state_abb, ~ fetch_gtrend(.)))

gtrend_state_flu_merge <- unnest(gtrend_state_list) %>%
  inner_join(mutate(state_matchup, state_abb = toupper(abb)), by = "state_abb") %>%
  select(-state_abb, -abb, location = name)

gtrend_merge <- bind_rows(gtrend_US_flu_merge, gtrend_state_flu_merge)

saveRDS(gtrend_merge, "Data/gtrend.Rds")
