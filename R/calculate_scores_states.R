rm(list = ls())
start_time <- Sys.time()

# Generate scores
library(tidyverse)
library(FluSight)
library(cdcfluview)

source("R/utils.R")
source('R/EpiDataAPI.R')

# Read in retrospective forecasts -------
forecasts_1415 <- read_forecasts("State Forecasts/Training/2014-2015", challenge = "state_ili")
forecasts_1516 <- read_forecasts("State Forecasts/Training/2015-2016", challenge = "state_ili")
forecasts_1617 <- read_forecasts("State Forecasts/Training/2016-2017", challenge = "state_ili")
forecasts_1718 <- read_forecasts("State Forecasts/Training/2017-2018", challenge = "state_ili")
forecasts_1819 <- read_forecasts("State Forecasts/Training/2018-2019", challenge = "state_ili")
forecasts_1920 <- read_forecasts("State Forecasts/Training/2019-2020", challenge = "state_ili")

# Create observed ILI -----
ili_init_pub_list <- readRDS('Data/ili_init_pub_list.Rds')

ILI_1415 <- ili_init_pub_list[["201740"]] %>%
  filter(week >= 40 | week < 25, season == "2014/2015",
         location %in% state.name) %>%
  mutate(week = week_inorder(week, season)) %>%
  arrange(location, week) %>%
  mutate(week = week_reset(week, season)) %>%
  select(location, week, ILI)

ILI_1516 <- ili_init_pub_list[["201740"]] %>%
  filter(week >= 40 | week < 25, season == "2015/2016",
         location %in% state.name) %>%
  mutate(week = week_inorder(week, season)) %>%
  arrange(location, week) %>%
  mutate(week = week_reset(week, season)) %>%
  select(location, week, ILI)

ILI_1617 <-  ili_init_pub_list[["201740"]] %>%
  filter(week >= 40 | week < 25, season == "2016/2017",
         location %in% state.name) %>%
  mutate(week = week_inorder(week, season)) %>%
  arrange(location, week) %>%
  mutate(week = week_reset(week, season)) %>%
  select(location, week, ILI)

ILI_1718 <-  ili_init_pub_list[["201828"]] %>%
  filter(week >= 40 | week < 25, season == "2017/2018",
         location %in% state.name) %>%
  mutate(week = week_inorder(week, season)) %>%
  arrange(location, week) %>%
  mutate(week = week_reset(week, season)) %>%
  select(location, week, ILI)

ILI_1819 <-  ili_init_pub_list[["201928"]] %>%
  filter(week >= 40 | week < 25, season == "2018/2019",
         location %in% state.name) %>%
  mutate(week = week_inorder(week, season)) %>%
  arrange(location, week) %>%
  mutate(week = week_reset(week, season)) %>%
  select(location, week, ILI)

ILI_1920 <-  ili_init_pub_list[["202028"]] %>%
  filter(week >= 40 | week < 25, season == "2019/2020",
         location %in% state.name) %>%
  mutate(week = week_inorder(week, season)) %>%
  arrange(location, week) %>%
  mutate(week = week_reset(week, season)) %>%
  select(location, week, ILI)

# Create truth ----------------------------------------------------------------
truth_1415 <- create_truth(fluview = FALSE, year = 2014, weekILI = ILI_1415,
                           challenge = "state_ili", start_wk = 40, end_wk = 20)

truth_1516 <- create_truth(fluview = FALSE, year = 2015, weekILI = ILI_1516,
                           challenge = "state_ili", start_wk = 40, end_wk = 20)

truth_1617 <- create_truth(fluview = FALSE, year = 2016, weekILI = ILI_1617,
                           challenge = "state_ili", start_wk = 40, end_wk = 20)

truth_1718 <- create_truth(fluview = FALSE, year = 2017, weekILI = ILI_1718,
                           challenge = "state_ili", start_wk = 40, end_wk = 20)

truth_1819 <- create_truth(fluview = FALSE, year = 2018, weekILI = ILI_1819,
                           challenge = "state_ili", start_wk = 40, end_wk = 20)

truth_1920 <- create_truth(fluview = FALSE, year = 2019, weekILI = ILI_1920,
                           challenge = "state_ili", start_wk = 40, end_wk = 20)

# # Expand observed truth to include all bins that will be counted as correct ----
# exp_truth_1415 <- expand_truth(truth_1415, week_expand = 1, percent_expand = 5,
#                                challenge = "state_ili")
# 
# exp_truth_1516 <- expand_truth(truth_1516, week_expand = 1, percent_expand = 5,
#                                challenge = "state_ili")
# 
# exp_truth_1617 <- expand_truth(truth_1617, week_expand = 1, percent_expand = 5,
#                                challenge = "state_ili")
# 
# exp_truth_1718 <- expand_truth(truth_1718, week_expand = 1, percent_expand = 5,
#                                challenge = "state_ili")
# 
# exp_truth_1819 <- expand_truth(truth_1819, week_expand = 1, percent_expand = 5,
#                                challenge = "state_ili")

# Score entries -----
full_scores_1415 <- calc_scores(forecasts_1415, truth_1415, 
                                season = "2014/2015", exclude = FALSE, 
                                eval = FALSE)

full_scores_1516 <- calc_scores(forecasts_1516, truth_1516, 
                                season = "2015/2016", exclude = FALSE, 
                                eval = FALSE)

full_scores_1617 <- calc_scores(forecasts_1617, truth_1617, 
                                season = "2016/2017", exclude = FALSE, 
                                eval = FALSE)

full_scores_1718 <- calc_scores(forecasts_1718, truth_1718, 
                                season = "2017/2018", exclude = FALSE, 
                                eval = FALSE)

full_scores_1819 <- calc_scores(forecasts_1819, truth_1819, 
                                season = "2018/2019", exclude = FALSE, 
                                eval = FALSE)

full_scores_1920 <- calc_scores(forecasts_1920, truth_1920, 
                                season = "2019/2020", exclude = FALSE, 
                                eval = FALSE)

state_full_scores <- bind_rows(full_scores_1415, full_scores_1516, 
                               full_scores_1617, full_scores_1718,
                               full_scores_1819, full_scores_1920)

# Save scores -----
saveRDS(state_full_scores, file = "Data/state_full_model_scores.Rds")


# Determine best fitting model by state
scores_by_state <- state_full_scores %>%
  filter(target %in% c("1 wk ahead", "2 wk ahead", "3 wk ahead", "4 wk ahead"),
         team != "ens-optimal-state") %>%
  group_by(team, location) %>%
  summarize(
    avg_score = mean(score),
    min_score = min(score)
  ) %>%
  group_by(location) %>%
  filter(avg_score == max(avg_score)) %>%
  ungroup() %>%
  select(team, location)

saveRDS(scores_by_state, file = "Data/state_best_model_fits.Rds")