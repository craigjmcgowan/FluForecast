# Generate scores
library(tidyverse)
library(FluSight)
library(cdcfluview)

source("R/utils.R")

# Read in retrospective forecasts -------
forecasts_1920 <- read_forecasts("CDC Submissions/2018-2019")

# Create observed ILI -----
ILI_1920 <- readRDS("Data/ili_init_pub_list.rds")[['202028']] %>%
  filter(!location %in% state.name) %>%
  mutate(order_week = week_inorder(week, "2019/2020")) %>%
  arrange(location, order_week) %>%
  select(location, week, ILI) %>%
  mutate(ILI = round(ILI, 1))

# Create truth ----------------------------------------------------------------
truth_1920 <- create_truth(fluview = FALSE, year = 2019, weekILI = ILI_1920,
                           challenge = "ilinet", start_wk = 42, end_wk = 18)

# Create evaluation period -----
eval_period_1920 <- create_eval_period(ILI_1920, truth_1920, "2019/2020")


# Score entries -----
full_scores_1920 <- calc_scores(forecasts_1920, truth_1920, 
                                season = "2019/2020", exclude = FALSE, 
                                eval = FALSE)

eval_scores_1920 <- calc_scores(forecasts_1920, truth_1920, 
                                season = "2019/2020", exclude = FALSE, 
                                eval = TRUE, eval_period = eval_period_1920)


# Save scores -----
saveRDS(eval_scores_1920, "Data/cdc_model_scores_1920.RDS")

# Look at different groupings of scores -----
group_by(full_scores_1920, team) %>%
  summarize(mean(score))

group_by(full_scores_1920, team, target) %>%
  summarize(mean(score))

filter(full_scores_1920, forecast_week < 40, forecast_week > 4) %>%
  group_by(team, target) %>%
  summarize(mean(score))
