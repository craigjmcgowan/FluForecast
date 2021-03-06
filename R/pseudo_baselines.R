# Create pseudo-baselines for calculating onset prior to 2007
# Calculating pseudo-baselines for 2003, 2004, 2005, and 2006 flu seasons
# Including 2007 as a test
library(tidyverse)

source("R/utils.R")

# Load data -------
ili_current <- readRDS("Data/ili_current.RDS")
virologic_combined <- readRDS("Data/virologic.RDS")


# Only keep years with full season data and before onset was officially calculated
ILI_collapse <- ili_current %>%
  filter(season %in% c("2002/2003", "2003/2004", "2004/2005", "2005/2006", "2006/2007"),
         !location %in% state.name) %>%
  # filter(season %in% c("2014/2015", "2015/2016", "2016/2017")) %>%
  select(location, season, week, year, ILI) %>%
  left_join(virologic_combined %>%
              filter(season %in% c("2002/2003", "2003/2004", "2004/2005", "2005/2006", "2006/2007")) %>%
              select(location, season, week, year, total_specimens,
                     percent_positive),
            by = c("location", "season", "week", "year")) %>%
  mutate(pos_specimens = round(total_specimens * percent_positive / 100),
         order_week = week_inorder(week, season)) %>%
  select(-percent_positive) %>%
  # Calculate total number of positive specimens in a season
  group_by(location, season) %>%
  mutate(tot_positive = sum(pos_specimens)) %>%
  ungroup() %>%
  # Calculate each week's percent of total positive specimens
  mutate(percent_total = pos_specimens / tot_positive) %>%
  # Flag weeks to be included in baseline calculations %>%
  group_by(location) %>%
  arrange(season, order_week, .by_group = TRUE) %>%
  mutate(lag_percent = lag(percent_total),
         lead_percent = lead(percent_total),
         # non_season_week = ifelse(
         #   (percent_total <= 0.02 & lag_percent <= 0.02) &
         #     (percent_total <= 0.02 & lead_percent <= 0.02),
         #   1, 0
         #  )
         non_season_week = ifelse(percent_total <= 0.02 & 
                                    (lag_percent <= 0.02 | lead_percent <= 0.02), 1, 0)
         ) %>%
  ungroup() %>%
  # Only keep non-season weeks
  filter(non_season_week == 1)

# # 2017 baseline as test - close enough with only keeping individual weeks
# analysis %>%
#   # filter(season %in% c("2004/2005", "2005/2006", "2006/2007")) %>%
#   group_by(location) %>%
#   summarize(mean_ili = mean(ILI),
#             sd_ili = sd(ILI),
#             calc_baseline = round(mean_ili + 2*sd_ili, 1)) %>%
#   left_join(FluSight::past_baselines %>%
#               filter(year == 2017),
#             by = c("location"))

# Calculate season specific baselines --------
# 2003 baselines
pseudo_baselines <- ILI_collapse %>%
  # 2003 baselines
  filter(season %in% c("2002/2003")) %>%
  group_by(location) %>%
  summarize(baseline = round(mean(ILI) + 2*sd(ILI), 1),
            year = 2003) %>%
  # 2004 baselines
  bind_rows(ILI_collapse %>%
              filter(season %in% c("2002/2003", "2003/2004")) %>%
              group_by(location) %>%
              summarize(baseline = round(mean(ILI) + 2*sd(ILI), 1),
                        year = 2004)) %>%
  # 2005 baselines
  bind_rows(ILI_collapse %>%
              filter(season %in% c("2002/2003", "2003/2004", "2004/2005")) %>%
              group_by(location) %>%
              summarize(baseline = round(mean(ILI) + 2*sd(ILI), 1),
                        year = 2005)) %>%
  # 2006 baselines
  bind_rows(ILI_collapse %>%
              filter(season %in% c("2003/2004", "2004/2005", "2005/2006")) %>%
              group_by(location) %>%
              summarize(baseline = round(mean(ILI) + 2*sd(ILI), 1),
                        year = 2006)) %>%
  # Add season baselines are applicable for
  mutate(season = paste0(year, "/", year + 1)) %>%
  select(-year)

# Merge baselines into ILI ------
pseudo_onsets <- ili_current %>%
  filter(season %in% c("2003/2004", "2004/2005", "2005/2006", "2006/2007"),
         location != "HHS Region 10") %>%
  left_join(pseudo_baselines,
            by = c("location", "season")) %>%
  group_by(location, season) %>%
  do(create_pseudo_onset(.)) %>%
  ungroup() 

# HHS Region 10 data no longer available as of 6/2019 - use previously calculated
region10_onsets <- tibble(
  season = c("2003/2004", "2004/2005", "2005/2006", "2006/2007"), 
  target = c("Season onset", "Season onset", "Season onset", "Season onset"), 
  location = c("HHS Region 10", "HHS Region 10","HHS Region 10", "HHS Region 10"),
  forecast_week = c(NA_integer_, NA_integer_, NA_integer_, NA_integer_), 
  bin_start_incl = c("44.0", "52.0", "51.0", "2.0")
)
  
pseudo_onsets <- bind_rows(pseudo_onsets, region10_onsets)  

# Save results
saveRDS(pseudo_onsets, "Data/pseudo_onsets.RDS")
