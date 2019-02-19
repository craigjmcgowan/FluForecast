rm(list = ls())

# Generate scores
library(tidyverse)
library(FluSight)
library(cdcfluview)

source("R/utils.R")
source('R/EpiDataAPI.R')

# Read in retrospective forecasts -------
forecasts_1415 <- read_forecasts("State Forecasts/2014-2015", challenge = "state_ili")
forecasts_1516 <- read_forecasts("State Forecasts/2015-2016", challenge = "state_ili")
forecasts_1617 <- read_forecasts("State Forecasts/2016-2017", challenge = "state_ili")
forecasts_1718 <- read_forecasts("State Forecasts/2017-2018", challenge = "state_ili")
forecasts_1819 <- read_forecasts("State Forecasts/2018-2019", challenge = "state_ili")

# Create observed ILI -----
ILI_1415 <- Epidata$fluview(regions = list("pa", "de", "md",
                                           "va", "wv"),
                            epiweeks = list(Epidata$range(201440, 201528)),
                            issue = 201740)$epidata %>%
  modify_depth(2, function(x) ifelse(is.null(x), NA, x)) %>%
  bind_rows() %>%
  mutate(week = as.integer(substr(epiweek, 5, 6)),
         location = case_when(
           region == "pa" ~ "Pennsylvania",
           region == "de" ~ "Delaware",
           region == "md" ~ "Maryland",
           region == "va" ~ "Virginia",
           region == "wv" ~ "West Virginia"
         ),
         ili = round(ili, 1)) %>%
  select(location, week, ILI = ili)

ILI_1516 <- Epidata$fluview(regions = list("pa", "de", "md",
                                           "va", "wv"),
                            epiweeks = list(Epidata$range(201540, 201628)),
                            issue = 201740)$epidata %>%
  modify_depth(2, function(x) ifelse(is.null(x), NA, x)) %>%
  bind_rows() %>%
  mutate(week = as.integer(substr(epiweek, 5, 6)),
         location = case_when(
           region == "pa" ~ "Pennsylvania",
           region == "de" ~ "Delaware",
           region == "md" ~ "Maryland",
           region == "va" ~ "Virginia",
           region == "wv" ~ "West Virginia"
         ),
         ili = round(ili, 1)) %>%
  select(location, week, ILI = ili)

ILI_1617 <-  Epidata$fluview(regions = list("pa", "de", "md",
                                            "va", "wv"),
                             epiweeks = list(Epidata$range(201640, 201728)),
                             issue = 201740)$epidata %>%
  modify_depth(2, function(x) ifelse(is.null(x), NA, x)) %>%
  bind_rows() %>%
  mutate(week = as.integer(substr(epiweek, 5, 6)),
         location = case_when(
           region == "pa" ~ "Pennsylvania",
           region == "de" ~ "Delaware",
           region == "md" ~ "Maryland",
           region == "va" ~ "Virginia",
           region == "wv" ~ "West Virginia"
         ),
         ili = round(ili, 1)) %>%
  select(location, week, ILI = ili)

ILI_1718 <-  Epidata$fluview(regions = list("pa", "de", "md",
                                            "va", "wv"),
                             epiweeks = list(Epidata$range(201740, 201828)),
                             issue = 201828)$epidata %>%
  modify_depth(2, function(x) ifelse(is.null(x), NA, x)) %>%
  bind_rows() %>%
  mutate(week = as.integer(substr(epiweek, 5, 6)),
         location = case_when(
           region == "pa" ~ "Pennsylvania",
           region == "de" ~ "Delaware",
           region == "md" ~ "Maryland",
           region == "va" ~ "Virginia",
           region == "wv" ~ "West Virginia"
         ),
         ili = round(ili, 1)) %>%
  select(location, week, ILI = ili)

ILI_1819 <-  Epidata$fluview(regions = list("pa", "de", "md",
                                            "va", "wv"),
                             epiweeks = list(Epidata$range(201840, 201906)),
                             issue = 201906)$epidata %>%
  modify_depth(2, function(x) ifelse(is.null(x), NA, x)) %>%
  bind_rows() %>%
  mutate(week = as.integer(substr(epiweek, 5, 6)),
         location = case_when(
           region == "pa" ~ "Pennsylvania",
           region == "de" ~ "Delaware",
           region == "md" ~ "Maryland",
           region == "va" ~ "Virginia",
           region == "wv" ~ "West Virginia"
         ),
         ili = round(ili, 1)) %>%
  select(location, week, ILI = ili)

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


# Expand observed truth to include all bins that will be counted as correct ----
exp_truth_1415 <- expand_truth(truth_1415, week_expand = 1, percent_expand = 5,
                               challenge = "state_ili")

exp_truth_1516 <- expand_truth(truth_1516, week_expand = 1, percent_expand = 5,
                               challenge = "state_ili")

exp_truth_1617 <- expand_truth(truth_1617, week_expand = 1, percent_expand = 5,
                               challenge = "state_ili")

exp_truth_1718 <- expand_truth(truth_1718, week_expand = 1, percent_expand = 5,
                               challenge = "state_ili")

exp_truth_1819 <- expand_truth(truth_1819, week_expand = 1, percent_expand = 5,
                               challenge = "state_ili")

# Score entries -----
full_scores_1415 <- calc_scores(forecasts_1415, exp_truth_1415, 
                                season = "2014/2015", exclude = FALSE, 
                                eval = FALSE)

full_scores_1516 <- calc_scores(forecasts_1516, exp_truth_1516, 
                                season = "2015/2016", exclude = FALSE, 
                                eval = FALSE)

full_scores_1617 <- calc_scores(forecasts_1617, exp_truth_1617, 
                                season = "2016/2017", exclude = FALSE, 
                                eval = FALSE)

full_scores_1718 <- calc_scores(forecasts_1718, exp_truth_1718, 
                                season = "2017/2018", exclude = FALSE, 
                                eval = FALSE)

full_scores_1819 <- calc_scores(forecasts_1819, exp_truth_1819, 
                                season = "2018/2019", exclude = FALSE, 
                                eval = FALSE)

state_full_scores <- bind_rows(full_scores_1415, full_scores_1516, 
                               full_scores_1617, full_scores_1718,
                               full_scores_1819)

# Save scores -----
save(state_full_scores,
     file = "Data/state_model_scores.Rdata")


# Determine best fitting model by state
scores_by_state <- state_full_scores %>%
  filter(season != "2018/2019") %>%
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

save(scores_by_state,
     file = "Data/state_best_model_fits.Rdata")
