# Create pseudo-baselines for calculating onset prior to 2007
# Calculating pseudo-baselines for 2003, 2004, 2005, and 2006 flu seasons
# Including 2007 as a test
library(tidyverse)

# Load data -------
load("Data/ili.Rdata")
load("Data/virologic.Rdata")
source("R/week_order_functions.R")

# Create pseudo-onsets from calculated baselines --------
create_pseudo_onset <- function(weekILI) {
  
  # Save maximum MMWR week in season being analyzed
  maxMMWR <- max(weekILI$week)
  
  # Add 52/53 to weeks in new year to keep weeks in order
  weekILI$week[weekILI$week < 40] <-
    as.integer(weekILI$week[weekILI$week < 40] + maxMMWR)
  
  # Check to see if 3 weeks above baseline have passed
  j <- 0  # Counter for weeks above peak
  for (i in head(weekILI$week, n = 1):tail(weekILI$week, n = 1)) {
    if (weekILI$ILI[weekILI$week == i] >= weekILI$baseline[weekILI$week == i]) {
      j <- j + 1
    } else {
      j <- 0
    }
    if (j == 3) {
      onset <- i - 2
      break
    }
    if (i == tail(weekILI$week, n = 1)) {
      onset <- "none"
    }
  }
  
  # If onset week > 52, reset to MMWR week
  if (is.numeric(onset) && onset > maxMMWR) {
    onset <- onset - maxMMWR
  }
  
  onset_truth <- data.frame(target = "Season onset",
                            location = weekILI$location[1],
                            forecast_week = as.integer(NA),
                            bin_start_incl = as.character(onset),
                            stringsAsFactors = FALSE) %>%
    mutate(bin_start_incl = trimws(replace(bin_start_incl,
                                           !is.na(bin_start_incl) & bin_start_incl != "none",
                                           format(round(as.numeric(
                                             bin_start_incl[!is.na(bin_start_incl) & bin_start_incl != "none"])
                                             , 1), nsmall = 1, trim = T))))
  
  
  return(onset_truth)
}  

# Only keep years with full season data and before onset was officially calculated
ILI_collapse <- ili_current %>%
  filter(season %in% c("2002/2003", "2003/2004", "2004/2005", "2005/2006", "2006/2007")) %>%
  # filter(season %in% c("2014/2015", "2015/2016", "2016/2017")) %>%
  select(location, season, week, year, ILI) %>%
  left_join(virologic_combined %>%
              filter(season %in% c("2002/2003", "2003/2004", "2004/2005", "2005/2006", "2006/2007")) %>%
              # filter(season %in% c("2014/2015", "2015/2016", "2016/2017")) %>%
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
         non_season_week = ifelse(percent_total <= 0.02, 1, 0)
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
  filter(season %in% c("2003/2004", "2004/2005", "2005/2006", "2006/2007")) %>%
  left_join(pseudo_baselines,
            by = c("location", "season")) %>%
  group_by(location, season) %>%
  do(create_pseudo_onset(.)) %>%
  ungroup() 

  
  
  
# Save results
save(pseudo_baselines, pseudo_onsets,
     file = "Data/pseudo_onsets_2003_2006.Rdata")
