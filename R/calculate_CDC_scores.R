# Generate scores
library(tidyverse)
library(FluSight)
library(cdcfluview)

source("R/utils.R")

# Read in retrospective forecasts -------
forecasts_1819 <- read_forecasts("CDC Submissions/2018-2019")

# Create observed ILI -----

ILI_1819 <- readRDS("Data/ili_init_pub_list.rds")[['201928']] %>%
  filter(!location %in% state.name) %>%
  select(location, week, ILI)

# Create truth ----------------------------------------------------------------
truth_1819 <- read_csv("CDC Scores/Targets_2018-2019.csv")


# Expand observed truth to include all bins that will be counted as correct ----
exp_truth_1819 <- expand_truth(truth_1819, week_expand = 1, percent_expand = 5,
                               challenge = "ilinet")

# Create evaluation period -----
eval_period_1819 <- create_eval_period(ILI_1819, truth_1819, "2018/2019")
eval_check <- read_csv("CDC Scores/Eval_Period_Bounds_2018-2019.csv") %>%
  mutate(end_week = case_when(end_week < 40 ~ end_week + 52,
                              TRUE ~ end_week))



# Score entries -----
full_scores_1819 <- calc_scores(forecasts_1819, exp_truth_1819, 
                                season = "2018/2019", exclude = FALSE, 
                                eval = FALSE)

eval_scores_1819 <- calc_scores(forecasts_1819, exp_truth_1819, 
                                season = "2018/2019", exclude = FALSE, 
                                eval = TRUE, eval_period = eval_period_1819)


# Save scores -----
saveRDS(eval_scores_1819, "Data/cdc_model_scores.RDS")


eval_scores_1819 %>%
  group_by(team) %>%
  summarize(mean_score = mean(score),
            skill = exp(mean_score))

