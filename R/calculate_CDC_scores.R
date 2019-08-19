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
  mutate(order_week = week_inorder(week, "2018/2019")) %>%
  arrange(location, order_week) %>%
  select(location, week, ILI) %>%
  mutate(ILI = round(ILI, 1))

# Create truth ----------------------------------------------------------------
truth_1819 <- create_truth(fluview = FALSE, year = 2018, weekILI = ILI_1819,
             challenge = "ilinet", start_wk = 42, end_wk = 18)

# Expand observed truth to include all bins that will be counted as correct ----
exp_truth_1819 <- expand_truth(truth_1819, week_expand = 1, percent_expand = 5,
                               challenge = "ilinet")

# Create evaluation period -----
eval_period_1819 <- create_eval_period(ILI_1819, truth_1819, "2018/2019")


# Score entries -----
full_scores_1819 <- calc_scores(forecasts_1819, exp_truth_1819, 
                                season = "2018/2019", exclude = FALSE, 
                                eval = FALSE)

eval_scores_1819 <- calc_scores(forecasts_1819, exp_truth_1819, 
                                season = "2018/2019", exclude = FALSE, 
                                eval = TRUE, eval_period = eval_period_1819)


# Save scores -----
saveRDS(eval_scores_1819, "Data/cdc_model_scores.RDS")

