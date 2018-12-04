# Dynamic harmonic regression model 

library(tidyverse)
library(forecast)
library(lubridate)
library(FluSight)
library(MMWRweek)
library(multidplyr)

# Load functions
source("R/utils.R")

# Load data
load("Data/ili.Rdata")
load("Data/virologic.Rdata")
load("Data/Gtrends.Rdata")
load("Data/state_gtrend.Rdata")

# Create clusters for use later in program
cluster <- create_cluster(cores = parallel::detectCores())

# Create truth for all seasons and combine datasets -------
load("Data/truth_and_data.Rdata")
# nested_truth <- ili_current %>%
#   filter(year >= 2010, season != "2009/2010") %>%
#   select(season, location, week, ILI) %>%
#   nest(-season) %>%
#   mutate(truth = map2(season, data,
#                       ~ create_truth(fluview = FALSE, year = substr(.x, 1, 4),
#                                      weekILI = .y)),
#          eval_period = pmap(list(data, truth, season),
#                             ~ create_eval_period(..1, ..2, ..3)),
#          exp_truth = map(truth,
#                          ~ expand_truth(.))) %>%
#   select(-data)
# 
# flu_data_merge <- select(ili_current, epiweek, ILI, year, week, season, location) %>%
#   inner_join(select(virologic_combined, location, season, year, week, cum_h1per,
#                    cum_h3per, cum_bper),
#             by = c("location", "season", "year", "week")) %>%
#   # Add Google Trends data by location
#   right_join(bind_rows(mutate(gtrend_US_flu_merge, location = "US National"),
#                       mutate(gtrend_MA_flu_merge, location = "HHS Region 1"),
#                       mutate(gtrend_NY_flu_merge, location = "HHS Region 2"),
#                       mutate(gtrend_PA_flu_merge, location = "HHS Region 3"),
#                       mutate(gtrend_FL_flu_merge, location = "HHS Region 4"),
#                       mutate(gtrend_IL_flu_merge, location = "HHS Region 5"),
#                       mutate(gtrend_TX_flu_merge, location = "HHS Region 6"),
#                       mutate(gtrend_MO_flu_merge, location = "HHS Region 7"),
#                       mutate(gtrend_CO_flu_merge, location = "HHS Region 8"),
#                       mutate(gtrend_CA_flu_merge, location = "HHS Region 9"),
#                       mutate(gtrend_WA_flu_merge, location = "HHS Region 10")) %>%
#               rename(region_hits = hits) %>%
#               select(-date),
#             by = c("location", "season", "year", "week")) %>%
#   right_join(gtrend_US_flu_merge,
#             by = c("season", "year", "week")) %>%
#   full_join(select(ili_backfill, -orig_ILI, -final_ILI),
#             by = c("location", "season", "year", "week")) %>%
#   # Remove 2008/2009 and 2009/2010 seasons due to pandemic activity
#   filter(!season %in% c("2008/2009", "2009/2010")) %>%
#   # Remove week 33 in 2014 so all seasons have 52 weeks - minimal activity
#   filter(!(year == 2014 & week == 33)) %>%
#   mutate(order_week = case_when(
#     week < 40 & season == "2014/2015" ~ week + 53,
#     week < 40 ~ week + 52,
#     TRUE ~ week
#   )) %>%
#   # Replace missing backfill data with Gaussian random draw from same season/week observed values
  # mutate(
  #   sim_backfill = map2_dbl(location, week, 
  #                           ~ rnorm(1, 
  #                              ili_backfill_avg$avg_backfill[
  #                                ili_backfill_avg$location == .x &
  #                                  ili_backfill_avg$week == .y],
  #                              ili_backfill_avg$sd_backfill[
  #                                ili_backfill_avg$location == .x &
  #                                  ili_backfill_avg$week == .y])),
  #   backfill = case_when(
  #     !is.na(backfill) ~ backfill,
  #     TRUE ~ sim_backfill
  #     )
  #   ) %>%
  # select(-sim_backfill)
# 
# save(nested_truth, flu_data_merge, file = "Data/truth_and_data.Rdata")



### Fit region-specific SARIMA model for each test season ###

season_arima_fits <- tibble()

for(this_season in c("2014/2015")) {#, "2015/2016", "2016/2017", "2017/2018")) {
  
  temp_data <- filter(ili_current, year <= as.numeric(substr(this_season, 1, 4)),
                      season != this_season) %>%
    # Remove week 53s and 2008/2009 and 2009/2010 seasons
    filter(!(season %in% c("1997/1998", "2003/2004", "2008/2009", "2014/2015") &
               week == 53),
           !(season %in% c("2008/2009", "2009/2010"))) %>%
    select(location, ILI) %>%
    nest(-location) %>%
    mutate(ILI = ifelse(ILI == 0, log(0.1), log(ILI)),
           ILI_ts = map(data,
                        ~ ts(.$ILI, frequency = 52, start = c(1999, 40))),
           group = rep(1:4, length.out = nrow(.)),
           season = this_season)
  
  temp_party <- temp_data %>%
      partition(group, cluster = cluster)

  temp_party %>%
    cluster_library(c("tidyverse", "forecast"))

  temp_fits <- temp_party %>%
    mutate(arima_fit = map(ILI_ts,
                           ~ auto.arima(.))) %>%
    collect() %>%
    as.tibble() %>%
    ungroup() %>%
    select(season, location, arima_fit)
  
  season_arima_fits <- bind_rows(
    season_arima_fits,
    temp_fits
  )
  
  rm(temp_data, temp_party, temp_fits)
}


saveRDS(season_arima_fits, file = "Data/seasonal_arima_fits.Rds")

season_arima_fits <- readRDS("Data/seasonal_arima_fits.Rds")


# Set up data for model fitting in parallel
sarima_model_data <- crossing(season = c("2014/2015"), # "2015/2016", "2016/2017", "2017/2018"),
                              week = c(43:71),
                              location = unique(flu_data_merge$location)) %>%
  filter(week < 71 | season == "2014/2015") %>%
  mutate(epiweek = case_when(
    season == "2014/2015" & week > 53 ~ 
      as.numeric(paste0(substr(season, 6, 9), 
                        str_pad(week - 53, 2, "left", "0"))),
    week > 52 ~ 
      as.numeric(paste0(substr(season, 6, 9), 
                        str_pad(week - 52, 2, "left", "0"))),
    TRUE ~ as.numeric(paste0(substr(season, 1, 4), str_pad(week, 2, "left", "0")))
  )) %>%
  left_join(season_arima_fits, by = c("season", "location")) %>%
  # Set up grouping for parallel
  mutate(group = rep(1:length(cluster), length.out = nrow(.)))


# Set up party_df and load necessary libraries and functions 
sarima_by_group <- sarima_model_data %>%
  partition(group, cluster = cluster)

sarima_by_group %>% 
  cluster_library(c("tidyverse", "forecast", "lubridate", "FluSight", "MMWRweek")) %>% 
  cluster_assign_value("flu_data_merge", flu_data_merge) %>%
  cluster_assign_value("ili_init_pub_list", ili_init_pub_list) %>%
  cluster_assign_value("fit_to_forecast", fit_to_forecast) %>%
  cluster_assign_value("sample_predictive_trajectories_arima", 
                       sample_predictive_trajectories_arima)

# Create forecasts for different transformations
sarima_forecasts <- sarima_by_group %>%
  mutate(
    pred_data = pmap(list(season, week, location, epiweek), 
                     ~ filter(flu_data_merge, year <= as.numeric(substr(..1, 6, 9)),
                              season != paste0(substr(..1, 6, 9), "/",
                                               as.numeric(substr(..1, 6, 9)) + 1),
                              season != ..1 | order_week %in% 40:..2,
                              location == ..3) %>%
                       left_join(select(ili_init_pub_list[[paste(..4)]],
                                        ILI, epiweek, location),
                                 by = c("epiweek", "location")) %>%
                       mutate(ILI = ifelse(is.na(ILI.y), ILI.x, ILI.y)) %>%
                       select(-ILI.x, -ILI.y) %>%
                       mutate(ILI = ts(ILI, frequency = 52, start = c(2006, 40)))),
    max_week = ifelse(season == "2014/2015", 53, 52),
    # Fit models
    pred_fit = pmap(list(pred_data, arima_fit),
                    ~ Arima(..1$ILI, model = ..2)),
    pred_results = pmap(
      list(pred_fit, pred_data, location, season, max_week),
      ~ fit_to_forecast(object = ..1,
                        pred_data = ..2,
                        location = ..3,
                        season = ..4,
                        max_week = ..5,
                        npaths = 250)
    )
  ) %>%
  collect() %>%
  as.tibble()

sarima_forecasts <- sarima_forecasts %>% ungroup() %>%
  select(season, location, week, pred_results)

saveRDS(sarima_forecasts, file = "Data/sarima_forecasts.Rds") 

# Normalize probabilities and score forecasts 
sarima_scores <- sarima_forecasts %>% 
  mutate(pred_results = map2(pred_results, location,
                             ~ mutate(.x, location = .y) %>%
                               normalize_probs())) %>%
  select(season, week, pred_results) %>%
  unnest() %>%
  nest(-season, -week) %>%
  left_join(nested_truth, by = "season") %>%
  mutate(scores = map2(data, exp_truth,
                       ~ score_entry(.x, .y)),
         eval_scores = pmap(list(scores, eval_period, season),
                            ~ create_eval_scores(..1, ..2, ..3))) %>%
  select(season, eval_scores) %>%
  unnest() 

saveRDS(sarima_scores, file = "Data/sarima_scores.Rds")