# Pull in and save datasets from EpiData and other sources 
library(tidyverse)
library(cdcfluview)
library(MMWRweek)
library(gtrendsR)
library(lubridate)
library(FluSight)

# Load functions
source("R/utils.R")
source('R/EpiDataAPI.R')

# Save current MMWR week in format Epidata wants
current_week <- as.numeric(paste0(MMWRweek(Sys.Date())[[1]], MMWRweek(Sys.Date())[[2]] ))

# Fetch wILI data from EpiData API -------
# Save ILINet values for previous 26 weeks for each week's publication available
ili_init_pub_list <- list()

pull_initpub_epidata(200740)

for(i in c(201040:201052, 201101:201152, 201201:201252, 201301:201338,
           201340:201352, 201401:201453, 201501:201552, 201601:201652, 
           201701:201752, 201801:201840)) {
  ili_init_pub_list[[paste(i)]] <- pull_initpub_epidata(i) %>%
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
    select(-ili) %>%
    # Make ILI always > 0
    mutate(ILI = case_when(near(ILI, 0) ~ 0.1,
                           TRUE ~ ILI))
}

# Original ILI data for each week
ili_orig <- bind_rows(ili_init_pub_list) %>%
  filter(lag == 0)
  
# ILI currently ------
ili_current <- bind_rows(
  pull_curr_epidata(199740, 200139),
  pull_curr_epidata(200140, 200439),
  pull_curr_epidata(200440, 200739),
  pull_curr_epidata(200740, 201039),
  pull_curr_epidata(201040, 201339),
  pull_curr_epidata(201340, 201639),
  pull_curr_epidata(201640, current_week)) %>%
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
  select(-ili) %>%
  # Make ILI always > 0
  mutate(ILI = case_when(near(ILI, 0) ~ 0.1,
                         TRUE ~ ILI))

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

save(ili_orig, ili_current, ili_init_pub_list, ili_backfill, ili_backfill_avg,
     file = "Data/ili.Rdata")

# Fetch virologic data from CDC -----
virologic_national <- who_nrevss(region = "national", years = c(1997:2018))
virologic_region <- who_nrevss(region = "hhs", years = c(1997:2018))

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
  # Create mock H1 and H3 percents if all A samples had been tested
  rowwise() %>%
  mutate(pos_samples = sum(a_h1, a_2009_h1n1, a_h3, a_subtyping_not_performed, 
           a_unable_to_subtype, b, bvic, byam, na.rm = TRUE),
         h1_per_samples = (sum(a_h1, a_2009_h1n1, na.rm = TRUE) / 
                             sum(a_h1, a_2009_h1n1, a_h3, na.rm = TRUE)) *
           sum(a_h1, a_2009_h1n1, a_h3, a_subtyping_not_performed, 
               a_unable_to_subtype, na.rm = T) / pos_samples,
         h3_per_samples = a_h3 / sum(a_h1, a_2009_h1n1, a_h3, na.rm = TRUE) *
           sum(a_h1, a_2009_h1n1, a_h3, a_subtyping_not_performed, 
               a_unable_to_subtype, na.rm = T) / pos_samples,
         b_per_samples = sum(b, bvic, byam, na.rm = T) / pos_samples) %>%
  ungroup() %>%
  mutate_at(vars(c("h1_per_samples", "h3_per_samples", "b_per_samples")),
            function(x) ifelse(is.na(x), 0, x)) %>%
  select(location = region, season, year, week, a_h1, a_2009_h1n1, a_h3, a_subtyping_not_performed, 
         a_unable_to_subtype, b, bvic, byam, h1_per_samples, h3_per_samples, b_per_samples)

save(virologic_combined,
     file = "Data/virologic.Rdata")

# Fetch Google Trends data -----
US_flu_0407 <- gtrends(keyword = "flu",
                       geo = "US",
                       time = paste(MMWRweek2Date(2004, 41), MMWRweek2Date(2007, 40)))$interest_over_time %>%
  mutate(hits = case_when(hits == "<1" ~ 1,
                              TRUE ~ as.numeric(hits)))

US_flu_0611 <- gtrends(keyword = "flu",
                       geo = "US",
                       time = paste(MMWRweek2Date(2006, 41), MMWRweek2Date(2011, 40)))$interest_over_time %>%
  mutate(hits = case_when(hits == "<1" ~ 1,
                              TRUE ~ as.numeric(hits)))

US_flu_1015 <- gtrends(keyword = "flu",
                      geo = "US",
                      time = paste(MMWRweek2Date(2010, 41), MMWRweek2Date(2015, 40)))$interest_over_time %>%
  mutate(hits = case_when(hits == "<1" ~ 1,
                              TRUE ~ as.numeric(hits)))

US_flu_1419 <- gtrends(keyword = "flu",
                       geo = "US")$interest_over_time %>%
  mutate(hits = case_when(hits == "<1" ~ 1,
                              TRUE ~ as.numeric(hits)))

# Set up ratios to normalize all Gtrends data to scale from 2010-2015
gratio_0607 <- US_flu_ratio(US_flu_0407, US_flu_0611)
gratio_1011 <- US_flu_ratio(US_flu_0611, US_flu_1015)
inv_gratio_1415 <- US_flu_ratio(US_flu_1419, US_flu_1015)

# Merge Gtrends data and rescale to 2010-2015 scale
gtrend_US_flu_merge <- filter(US_flu_0407, !date %in% US_flu_0611$date) %>%
  mutate(hits = hits / gratio_0607) %>%
  bind_rows(US_flu_0611) %>%
  filter(!date %in% US_flu_1015$date) %>%
  mutate(hits = hits / gratio_1011) %>%
  bind_rows(US_flu_1015) %>%
  filter(!date %in% US_flu_1419$date) %>%
  bind_rows(US_flu_1419 %>%
              mutate(hits = hits / inv_gratio_1415)) %>%
  select(date, hits) %>%
  mutate(week = MMWRweek(date)[[2]],
         year = MMWRweek(date)[[1]],
         season = ifelse(week >= 40,
                         paste0(year, "/", year + 1),
                         paste0(year - 1, "/", year)),
         hits = case_when(near(hits, 0) ~ 1,
                          TRUE ~ hits))


save(gtrend_US_flu_merge, file = "Data/Gtrends.Rdata")
