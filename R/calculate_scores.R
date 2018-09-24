# Generate scores
library(tidyverse)
library(FluSight)
library(cdcfluview)

source("R/utils.R")

# Read in retrospective forecasts -------
forecasts_1011 <- read_forecasts("Forecasts/2010-2011")
forecasts_1112 <- read_forecasts("Forecasts/2011-2012")
forecasts_1213 <- read_forecasts("Forecasts/2012-2013")
forecasts_1314 <- read_forecasts("Forecasts/2013-2014")
forecasts_1415 <- read_forecasts("Forecasts/2014-2015")
forecasts_1516 <- read_forecasts("Forecasts/2015-2016")
forecasts_1617 <- read_forecasts("Forecasts/2016-2017")
forecasts_1718 <- read_forecasts("Forecasts/2017-2018")

# Create observed ILI -----
ILI_1011 <- ilinet(region = "national", year = 2010) %>%
  mutate(location = "US National") %>%
  select(location, week, ILI = weighted_ili) %>%
  bind_rows(
    ilinet(region = "hhs", year = 2010) %>%
      mutate(location = paste("HHS", region)) %>%
      select(location, week, ILI = weighted_ili)
  ) %>%
  mutate(ILI = round(ILI, 1))

ILI_1112 <- ilinet(region = "national", year = 2011) %>%
  mutate(location = "US National") %>%
  select(location, week, ILI = weighted_ili) %>%
  bind_rows(
    ilinet(region = "hhs", year = 2011) %>%
      mutate(location = paste("HHS", region)) %>%
      select(location, week, ILI = weighted_ili)
  ) %>%
  mutate(ILI = round(ILI, 1))

ILI_1213 <- ilinet(region = "national", year = 2012) %>%
    mutate(location = "US National") %>%
    select(location, week, ILI = weighted_ili) %>%
    bind_rows(
      ilinet(region = "hhs", year = 2012) %>%
        mutate(location = paste("HHS", region)) %>%
        select(location, week, ILI = weighted_ili)
    ) %>%
    mutate(ILI = round(ILI, 1))

ILI_1314 <- ilinet(region = "national", year = 2013) %>%
  mutate(location = "US National") %>%
  select(location, week, ILI = weighted_ili) %>%
  bind_rows(
    ilinet(region = "hhs", year = 2013) %>%
      mutate(location = paste("HHS", region)) %>%
      select(location, week, ILI = weighted_ili)
  ) %>%
  mutate(ILI = round(ILI, 1))

ILI_1415 <- ilinet(region = "national", year = 2014) %>%
  mutate(location = "US National") %>%
  select(location, week, ILI = weighted_ili) %>%
  bind_rows(
    ilinet(region = "hhs", year = 2014) %>%
      mutate(location = paste("HHS", region)) %>%
      select(location, week, ILI = weighted_ili)
  ) %>%
  mutate(ILI = round(ILI, 1))

ILI_1516 <- read_csv("Data/ILINet_US_wk28_2016.csv") %>%
  select(week = WEEK, location = REGION.TYPE, ILI = X..WEIGHTED.ILI) %>%
  filter(week >= 42 | week <= 22) %>%
  mutate(location = paste0("US ", location), ILI = round(ILI, 1)) %>%
  bind_rows(
    read_csv("Data/ILINet_Regional_wk28_2016.csv") %>%
      select(location = REGION, week = WEEK,  ILI = X..WEIGHTED.ILI) %>%
      filter(week >= 42 | week <= 22) %>%
      mutate(location = paste0("HHS ", location), ILI = round(ILI, 1))
  )

ILI_1617 <- read_csv(paste0("Data/ILINet_Regional_wk28_2017.csv")) %>%
  select(REGION, `% WEIGHTED ILI`, WEEK) %>%
  mutate(REGION = paste("HHS", REGION)) %>%
  bind_rows(read_csv(paste0("Data/ILINet_US_wk28_2017.csv")) %>%
              select(REGION, `% WEIGHTED ILI`, WEEK) %>%
              mutate(REGION = "US National")) %>%
  rename(location = REGION,
         ILI = `% WEIGHTED ILI`,
         week = WEEK) %>%
  mutate(ILI = round(ILI, 1))

ILI_1718 <- read_csv("Data/ILINet_US_wk28_2018.csv") %>%
  mutate(location = "US National") %>%
  select(location, week, ILI = weighted_ili) %>%
  bind_rows(
    read_csv("Data/ILINet_Regional_wk28_2018.csv") %>%
      mutate(location = paste("HHS", region)) %>%
      select(location, week, ILI = weighted_ili)
  ) %>%
  mutate(ILI = round(ILI, 1))

# Create truth ----------------------------------------------------------------
truth_1011 <- create_truth(fluview = FALSE, year = 2010, weekILI = ILI_1011,
                           challenge = "ilinet")

truth_1112 <- create_truth(fluview = FALSE, year = 2011, weekILI = ILI_1112,
                           challenge = "ilinet")

truth_1213 <- create_truth(fluview = FALSE, year = 2012, weekILI = ILI_1213,
                           challenge = "ilinet")

truth_1314 <- create_truth(fluview = FALSE, year = 2013, weekILI = ILI_1314,
                           challenge = "ilinet")

truth_1415 <- create_truth(fluview = FALSE, year = 2014, weekILI = ILI_1415,
                           challenge = "ilinet")

truth_1516 <- create_truth(fluview = FALSE, year = 2015, weekILI = ILI_1213,
                           challenge = "ilinet")

truth_1617 <- create_truth(fluview = FALSE, year = 2016, weekILI = ILI_1213,
                           challenge = "ilinet")

truth_1718 <- create_truth(fluview = FALSE, year = 2017, weekILI = ILI_1718,
                           challenge = "ilinet")


# Expand observed truth to include all bins that will be counted as correct ----
exp_truth_1011 <- expand_truth(truth_1011, week_expand = 1, percent_expand = 5,
                               challenge = "ilinet")

exp_truth_1112 <- expand_truth(truth_1112, week_expand = 1, percent_expand = 5,
                               challenge = "ilinet")

exp_truth_1213 <- expand_truth(truth_1213, week_expand = 1, percent_expand = 5,
                               challenge = "ilinet")

exp_truth_1314 <- expand_truth(truth_1314, week_expand = 1, percent_expand = 5,
                               challenge = "ilinet")

exp_truth_1415 <- expand_truth(truth_1415, week_expand = 1, percent_expand = 5,
                               challenge = "ilinet")

exp_truth_1516 <- expand_truth(truth_1516, week_expand = 1, percent_expand = 5,
                               challenge = "ilinet")

exp_truth_1617 <- expand_truth(truth_1617, week_expand = 1, percent_expand = 5,
                               challenge = "ilinet")

exp_truth_1718 <- expand_truth(truth_1718, week_expand = 1, percent_expand = 5,
                               challenge = "ilinet")

# Create evaluation period -----
eval_period_1011 <- create_eval_period(ILI_1011, truth_1011, "2010/2011")

eval_period_1112 <- create_eval_period(ILI_1112, truth_1112, "2011/2012")

eval_period_1213 <- create_eval_period(ILI_1213, truth_1213, "2012/2013")

eval_period_1314 <- create_eval_period(ILI_1314, truth_1314, "2013/2014")

eval_period_1415 <- create_eval_period(ILI_1415, truth_1415, "2014/2015")

eval_period_1516 <- create_eval_period(ILI_1516, truth_1516, "2015/2016")

eval_period_1617 <- create_eval_period(ILI_1617, truth_1617, "2016/2017")

eval_period_1718 <- create_eval_period(ILI_1718, truth_1718, "2017/2018")


# Score entries -----
source("R/calc_scores_helper.R")
eval_scores_1011 <- calc_scores(forecasts_1011, exp_truth_1011, 
                                season = "2010/2011", exclude = FALSE, 
                                eval = TRUE, eval_period = eval_period_1011)

eval_scores_1112 <- calc_scores(forecasts_1112, exp_truth_1112, 
                                season = "2011/2012", exclude = FALSE, 
                                eval = TRUE, eval_period = eval_period_1112)

eval_scores_1213 <- calc_scores(forecasts_1213, exp_truth_1213, 
                                season = "2012/2013", exclude = FALSE, 
                                eval = TRUE, eval_period = eval_period_1213)

eval_scores_1314 <- calc_scores(forecasts_1314, exp_truth_1314, 
                                season = "2013/2014", exclude = FALSE, 
                                eval = TRUE, eval_period = eval_period_1314)

eval_scores_1415 <- calc_scores(forecasts_1415, exp_truth_1415, 
                                season = "2014/2015", exclude = FALSE, 
                                eval = TRUE, eval_period = eval_period_1415)

eval_scores_1516 <- calc_scores(forecasts_1516, exp_truth_1516, 
                                season = "2015/2016", exclude = FALSE, 
                                eval = TRUE, eval_period = eval_period_1516)

eval_scores_1617 <- calc_scores(forecasts_1617, exp_truth_1617, 
                                season = "2016/2017", exclude = FALSE, 
                                eval = TRUE, eval_period = eval_period_1617)

eval_scores_1718 <- calc_scores(forecasts_1718, exp_truth_1718, 
                                season = "2017/2018", exclude = FALSE, 
                                eval = TRUE, eval_period = eval_period_1718)

all_eval_scores <- bind_rows(eval_scores_1011, eval_scores_1112, 
                             eval_scores_1213, eval_scores_1314,
                             eval_scores_1415, eval_scores_1516, 
                             eval_scores_1617, eval_scores_1718)

# Save scores -----
save(eval_scores_1011, eval_scores_1112, eval_scores_1213, eval_scores_1314,
     eval_scores_1415, eval_scores_1516, eval_scores_1617, eval_scores_1718,
     all_eval_scores,
     file = "model_scores.Rdata")

