start_time <- Sys.time()

# Pull in and save datasets from EpiData and other sources 
library(tidyverse)
library(cdcfluview)
library(MMWRweek)
library(gtrendsR)
library(lubridate)
library(FluSight)
library(epiforecast)
library(zoo)


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
                   "hhs7", "hhs8", "hhs9", "hhs10", 'dc', 'pr'),
           name = c("US National", "HHS Region 1", "HHS Region 2", "HHS Region 3",
                    "HHS Region 4", "HHS Region 5", "HHS Region 6",  "HHS Region 7", 
                    "HHS Region 8", "HHS Region 9", "HHS Region 10", 'District of Columbia',
                    'Puerto Rico'))
  )

region_matchup <- tibble(
  state = c(state.name, 'Puerto Rico', 'District of Columbia'),
  region = c('HHS Region 4', 'HHS Region 10', 'HHS Region 9', 'HHS Region 6', 'HHS Region 9',
             'HHS Region 8', 'HHS Region 1', 'HHS Region 3', 'HHS Region 4', 'HHS Region 4',
             'HHS Region 9', 'HHS Region 10', 'HHS Region 5', 'HHS Region 5', 'HHS Region 7',
             'HHS Region 7', 'HHS Region 4', 'HHS Region 6', 'HHS Region 1', 'HHS Region 3', 
             'HHS Region 1', 'HHS Region 5', 'HHS Region 5', 'HHS Region 4', 'HHS Region 7',
             'HHS Region 8', 'HHS Region 7', 'HHS Region 9', 'HHS Region 1', 'HHS Region 2',
             'HHS Region 6', 'HHS Region 2', 'HHS Region 4', 'HHS Region 8', 'HHS Region 5',
             'HHS Region 6', 'HHS Region 10', 'HHS Region 3', 'HHS Region 1', 'HHS Region 4',
             'HHS Region 8', 'HHS Region 4', 'HHS Region 6', 'HHS Region 8', 'HHS Region 1',
             'HHS Region 3', 'HHS Region 10', 'HHS Region 3', 'HHS Region 5', 'HHS Region 8',
             'HHS Region 2', 'HHS Region 3'))

# ILI data ------
current_epidata <- function(start_wk, end_wk) {
  Epidata$fluview(
    regions = list("al", "ak", "az", "ar", "ca", "co", "ct", "de", "ga", "hi", 
                   "id", "il", "in", "ia", "ks", "ky", "la", "me", "md", "ma", 
                   "mi", "mn", "ms", "mo", "mt", "ne", "nv", "nh", "nj", "nm", 
                   "ny", "nc", "nd", "oh", "ok", "or", "pa", "ri", "sc", "sd", 
                   "tn", "tx", "ut", "vt", "va", "wa", "wv", "wi", "wy", "nat",
                   "hhs1", "hhs2", "hhs3", "hhs4", "hhs5", "hhs6", "hhs7", "hhs8",
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
  current_epidata(201840, 201939),
  current_epidata(201940, pull_week)
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
  summarize(avg_backfill = mean(backfill),
            sd_backfill = sd(backfill))

ili_backfill_month_avg <- ili_backfill %>%
  mutate(month = month(MMWRweek2Date(year, week))) %>%
  group_by(location, month) %>%
  summarize(avg_backfill = mean(backfill),
            sd_backfill = sd(backfill))

saveRDS(ili_orig, "Data/ili_orig.RDS")
saveRDS(ili_current, "Data/ili_current.RDS")
saveRDS(ili_init_pub_list, "Data/ili_init_pub_list.RDS")
saveRDS(ili_backfill, "Data/ili_backfill.RDS")
saveRDS(ili_backfill_avg, "Data/ili_backfill_avg.RDS")
saveRDS(ili_backfill_month_avg, "Data/ili_backfill_month_avg.RDS")

# Fetch virologic data from CDC -----
virologic_national <- who_nrevss(region = "national", years = c(1997:2020))
virologic_region <- who_nrevss(region = "hhs", years = c(1997:2020))
virologic_state <- who_nrevss(region = "state", years = c(2010:2020))

virologic_before_1516 <- bind_rows(
  virologic_national[[1]],
  virologic_region[[1]], 
  virologic_state[[1]] %>%
    mutate_at(c("a_2009_h1n1", "a_h1", "a_h3", "a_subtyping_not_performed",
                "a_unable_to_subtype", "b", "h3n2v", "total_specimens",
                "percent_positive"),
              as.integer) %>%
    filter(region %in% state.name)
)

virologic_ph_lab <- bind_rows(
  virologic_national[[2]], 
  virologic_region[[2]]
) 
 
# virologic_ph_lab_state <- virologic_state[[2]] %>%
#   mutate(season = paste(substr(season_description, 8, 11), 
#                         as.numeric(substr(season_description, 8, 11)) + 1,
#                         sep = "/")) %>%
#   mutate_at(c("a_2009_h1n1", "a_h3", "a_subtyping_not_performed",
#               "b", "bvic", "byam", "h3n2v", "total_specimens"),
#             as.integer) %>%
#   select(-region_type, -season_description, -wk_date) %>%
#   mutate_at(vars(c("a_2009_h1n1", "a_h3", "a_subtyping_not_performed",
#                    "b", "bvic", "byam", "h3n2v")),
#             function(x) ifelse(is.na(x), 0, x)) %>%
#   # Create mock H1 and H3 percents if all A samples had been tested
#   rowwise() %>%
#   mutate(h1sum = sum(a_2009_h1n1),
#          h3sum = sum(a_h3, h3n2v),
#          asum = sum(a_2009_h1n1, a_h3, a_subtyping_not_performed, h3n2v),
#          bsum = sum(b, bvic, byam),
#          # Subtype percentages of typed As
#          h1per_of_a = ifelse(is.na(h1sum / (h1sum + h3sum)), 0.5, 
#                              h1sum / (h1sum + h3sum)),
#          h3per_of_a = ifelse(is.na(h3sum / (h1sum + h3sum)), 0.5, 
#                              h3sum / (h1sum + h3sum)),
#          # Cumulative percentages of each type, assuming A type %s are representative
#          cum_bper = bsum / (asum + bsum),
#          cum_h1per = h1per_of_a * asum / (asum + bsum),
#          cum_h3per = h3per_of_a * asum / (asum + bsum),
#          # If no samples reported, make each type 33%
#          cum_h1per = ifelse(is.na(cum_h1per), 1/3, cum_h1per),
#          cum_h3per = ifelse(is.na(cum_h3per), 1/3, cum_h3per),
#          cum_bper = ifelse(is.na(cum_bper), 1/3, cum_bper)) %>% 
#   ungroup() %>%
#   # Create dataset with all weeks, with each week having the cumulative measure
#   full_join(crossing(season = c('2015/2016', '2016/2017', '2017/2018', '2018/2019'),
#                      week = c(1:52)),
#             by = "season") %>%
#   # Add year variable for future use
#   mutate(year = case_when(week < 40 ~ as.numeric(substr(season, 6, 9)),
#                           TRUE ~ as.numeric(substr(season, 1, 4)))) %>%
#   # Select variables of interest
#   select(location = region, season, week, year, a_2009_h1n1, a_h3, 
#          a_subtyping_not_performed, b, bvic, byam, 
#          h1per_of_a, h3per_of_a, cum_h1per, cum_h3per, cum_bper, total_specimens)
#   
# saveRDS(virologic_ph_lab_state, "Data/virologic_ph_lab_state.RDS")  


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
  # Create cumulative and rolling week influenza percentage measures
  group_by(season, region) %>%
  arrange(wk_date, .by_group = TRUE) %>%
  mutate(h1sum_6wk = rollapply(a_2009_h1n1, 6, sum, align = 'right', partial = TRUE) + 
           rollapply(a_h1, 6, sum, align = 'right', partial = TRUE),
         h3sum_6wk = rollapply(a_h3, 6, sum, align = 'right', partial = TRUE),
         asum_6wk = rollapply(a_pos_samples, 6, sum, align = 'right', partial = TRUE),
         bsum_6wk = rollapply(b_pos_samples, 6, sum, align = 'right', partial = TRUE),
         h1sum = cumsum(a_2009_h1n1) + cumsum(a_h1),
         h3sum = cumsum(a_h3),
         asum = cumsum(a_pos_samples),
         bsum = cumsum(b_pos_samples),
         # Subtype percentages of typed As
         h1per_of_a_6wk = ifelse(is.na(h1sum_6wk / (h1sum_6wk + h3sum_6wk)), 0.5, 
                                 h1sum_6wk / (h1sum_6wk + h3sum_6wk)),
         h3per_of_a_6wk = ifelse(is.na(h3sum_6wk / (h1sum_6wk + h3sum_6wk)), 0.5, 
                                 h3sum_6wk / (h1sum_6wk + h3sum_6wk)),
         h1per_of_a = ifelse(is.na(h1sum / (h1sum + h3sum)), 0.5, 
                             h1sum / (h1sum + h3sum)),
         h3per_of_a = ifelse(is.na(h3sum / (h1sum + h3sum)), 0.5, 
                             h3sum / (h1sum + h3sum)),
         # Cumulative percentages of each type, assuming A type %s are representative
         cum_bper_6wk = bsum_6wk / (asum_6wk + bsum_6wk),
         cum_h1per_6wk = h1per_of_a_6wk * asum_6wk / (asum_6wk + bsum_6wk),
         cum_h3per_6wk = h3per_of_a_6wk * asum_6wk / (asum_6wk + bsum_6wk),
         cum_bper = bsum / (asum + bsum),
         cum_h1per = h1per_of_a * asum / (asum + bsum),
         cum_h3per = h3per_of_a * asum / (asum + bsum)) %>%
  # If no samples reported, make each type 33%
  mutate(cum_h1per_6wk = ifelse(is.na(cum_h1per_6wk), 1/3, cum_h1per_6wk),
         cum_h3per_6wk = ifelse(is.na(cum_h3per_6wk), 1/3, cum_h3per_6wk),
         cum_bper_6wk = ifelse(is.na(cum_bper_6wk), 1/3, cum_bper_6wk),
         cum_h1per = ifelse(is.na(cum_h1per), 1/3, cum_h1per),
         cum_h3per = ifelse(is.na(cum_h3per), 1/3, cum_h3per),
         cum_bper = ifelse(is.na(cum_bper), 1/3, cum_bper)) %>% 
  ungroup() %>%
  select(location = region, season, year, week, a_h1, a_2009_h1n1, a_h3, 
         a_subtyping_not_performed, a_unable_to_subtype, b, bvic, byam, 
         h1per_of_a_6wk, h1per_of_a, h3per_of_a_6wk, h3per_of_a,
         cum_h1per_6wk, cum_h1per, cum_h3per_6wk, cum_h3per, 
         cum_bper_6wk, cum_bper, total_specimens, percent_positive)

# For 2015/2016 onwards, assign states the values from their HHS Region
# Data aren't published at the state level weekly :-/
virologic_state_ph_lab <- filter(virologic_combined, location != "US National",
                                 year >= 2015, season != '2014/2015') %>%
  full_join(region_matchup, by = c('location' = 'region')) %>%
  select(-location) %>%
  rename(location = state)

saveRDS(bind_rows(virologic_combined, virologic_state_ph_lab), file = "Data/virologic.Rds")

# Fetch Google Trends data -----
gtrend_US_flu_merge <- fetch_gtrend("US") %>%
  mutate(location = "US National")

gtrend_state_list <- tibble(state_abb = state.abb) %>%
  mutate(data = map(state_abb, ~ fetch_gtrend(.)))

gtrend_state_flu_merge <- unnest(gtrend_state_list, col = c(data)) %>%
  inner_join(mutate(state_matchup, state_abb = toupper(abb)), by = "state_abb") %>%
  select(-state_abb, -abb, location = name)

gtrend <- bind_rows(gtrend_US_flu_merge, gtrend_state_flu_merge)

saveRDS(gtrend, "Data/gtrend.Rds")

# COVID data -----
### Covidtracking.com
covid_tracking_us <- read_csv('https://covidtracking.com/api/v1/us/daily.csv') %>%
  mutate(location = 'US National',
         date = ymd(date),
         week = MMWRweek(date)$MMWRweek,
         year = year(date)) %>%
  group_by(week, year, location) %>%
  arrange(date, .by_group = TRUE) %>%
  summarize(date = first(date),
            positive = coalesce(last(positive), 0),
            negative = coalesce(last(negative), 0),
            recovered = coalesce(last(recovered), 0),
            death = coalesce(last(death), 0),
            positiveIncrease = coalesce(mean(positiveIncrease, na.rm = T), 0),
            negativeIncrease = coalesce(mean(negativeIncrease, na.rm = T), 0)) %>%
  ungroup()

covid_tracking_state_raw <- read_csv('https://covidtracking.com/api/v1/states/daily.csv') %>%
  rename(state_abb = state) %>%
  inner_join(mutate(state_matchup, state_abb = toupper(abb)), by = "state_abb") %>%
  inner_join(region_matchup, by = c('name' = 'state'))

covid_tracking_hhs <- covid_tracking_state_raw %>%
  group_by(date, region) %>%
  summarize(positive = coalesce(sum(positive, na.rm = T), 0),
            negative = coalesce(sum(negative, na.rm = T), 0),
            recovered = coalesce(sum(recovered, na.rm = T), 0),
            death = coalesce(sum(death, na.rm = T), 0),
            positiveIncrease = coalesce(sum(positiveIncrease, na.rm = T), 0),
            negativeIncrease = coalesce(sum(negativeIncrease, na.rm = T), 0)) %>%
  mutate(date = ymd(date),
         week = MMWRweek(date)$MMWRweek,
         year = year(date),
         location = region) %>%
  group_by(week, year, location) %>%
  arrange(date, .by_group = TRUE) %>%
  summarize(date = first(date),
            positive = coalesce(last(positive), 0),
            negative = coalesce(last(negative), 0),
            recovered = coalesce(last(recovered), 0),
            death = coalesce(last(death), 0),
            positiveIncrease = coalesce(mean(positiveIncrease, na.rm = T), 0),
            negativeIncrease = coalesce(mean(negativeIncrease, na.rm = T), 0)) %>%
  ungroup()

covid_tracking_state <- covid_tracking_state_raw %>%
  mutate(date = ymd(date),
         week = MMWRweek(date)$MMWRweek,
         year = year(date),
         location = name) %>%
  group_by(week, year, location) %>%
  arrange(date, .by_group = TRUE) %>%
  summarize(date = first(date),
            positive = coalesce(last(positive), 0),
            negative = coalesce(last(negative), 0),
            recovered = coalesce(last(recovered), 0),
            death = coalesce(last(death), 0),
            positiveIncrease = coalesce(mean(positiveIncrease, na.rm = T), 0),
            negativeIncrease = coalesce(mean(negativeIncrease, na.rm = T), 0)) %>%
  ungroup()
  
covid_tracking <- bind_rows(
  covid_tracking_us,
  covid_tracking_hhs,
  covid_tracking_state
)

saveRDS(covid_tracking, "Data/covid_tracking.RDS")

# Gtrends
covid_gtrend_us <- full_join(
  fetch_gtrend('US', 'covid') %>%
    rename('covid_hits' = 'hits'),
  fetch_gtrend('US', 'coronavirus') %>%
    select(date, coronavirus_hits = hits),
  by = 'date'
) %>%
  mutate(location = 'US National')

covid_gtrend_state <- tibble(state_abb = c('CA', 'WA', 'NY')) %>%
  mutate(data = map(state_abb, fetch_gtrend, 'covid')) %>%
  unnest(col = c(data)) %>%
  inner_join(mutate(state_matchup, state_abb = toupper(abb)), by = "state_abb") %>%
  select(-state_abb, -abb, location = name)

  
gtrend_state_list <- tibble(state_abb = state.abb) %>%
  mutate(data = map(state_abb, ~ fetch_gtrend(., 'coronavirus')))


saveRDS(covid_gtrend_us, "Data/covid_gtrend.RDS")

Sys.time() - start_time
