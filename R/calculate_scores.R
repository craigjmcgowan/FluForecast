# Generate scores
library(tidyverse)
library(FluSight)
library(cdcfluview)

source("R/read_forecasts.R")

# Read in retrospective forecasts -------
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
truth_1213 <- create_truth(fluview = FALSE, year = 2012, weekILI = ILI_1213,
                           challenge = "ilinet")

truth_1718 <- create_truth(fluview = FALSE, year = 2017, weekILI = ILI_1718,
                           challenge = "ilinet")

# Expand observed truth to include all bins that will be counted as correct
exp_truth_1213 <- expand_truth(truth_1213, week_expand = 1, percent_expand = 5,
                               challenge = "ilinet")

exp_truth_1718 <- expand_truth(truth_1718, week_expand = 1, percent_expand = 5,
                               challenge = "ilinet")
# Create evaluation period
eval_period_1213 <- create_eval_period(ILI_1213, truth_1213, "2012/2013")

eval_period_1718 <- create_eval_period(ILI_1718, truth_1718, "2017/2018")

# Score entries -----
source("R/calc_scores_helper.R")

eval_scores_1213 <- calc_scores(forecasts_1213, exp_truth_1213, 
                                season = "2012/2013", exclude = FALSE, 
                                eval = TRUE, eval_period = eval_period_1213)

eval_scores_1718 <- calc_scores(forecasts_1718, exp_truth_1718, 
                                season = "2017/2018", exclude = FALSE, 
                                eval = TRUE, eval_period = eval_period_1718)

# Save scores
save(eval_scores_1213, eval_scores_1718,
     file = "model_scores.Rdata")


overall_avg <- function(scores, this_location = NULL) {
  if (is.null(this_location)) this_location = unique(scores$location)
  scores %>%
    filter(location %in% this_location) %>%
    group_by(team) %>%
    summarize(score = mean(score, na.rm = T)) %>%
    mutate(skill = exp(score)) %>%
    arrange(desc(score))
}

overall_avg(eval_scores_1718)
