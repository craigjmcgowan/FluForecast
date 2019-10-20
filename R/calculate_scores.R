# Generate scores
library(tidyverse)
library(FluSight)
library(cdcfluview)

source("R/utils.R")

# Read in retrospective forecasts -------
forecasts_1011 <- read_forecasts("Forecasts/Training/2010-2011")
forecasts_1112 <- read_forecasts("Forecasts/Training/2011-2012")
forecasts_1213 <- read_forecasts("Forecasts/Training/2012-2013")
forecasts_1314 <- read_forecasts("Forecasts/Training/2013-2014")
forecasts_1415 <- read_forecasts("Forecasts/Training/2014-2015")
forecasts_1516 <- read_forecasts("Forecasts/Training/2015-2016")
forecasts_1617 <- read_forecasts("Forecasts/Training/2016-2017")
forecasts_1718 <- read_forecasts("Forecasts/Training/2017-2018")
forecasts_1819 <- read_forecasts("Forecasts/Training/2018-2019")

# Create observed ILI -----
ili_current <- readRDS('Data/ili_current.RDS') %>%
  filter(!location %in% state.name)

ILI_1011 <- ili_current %>%
  filter(season == "2010/2011") %>%
  select(season, location, week, ILI) %>%
  mutate(ILI = round(ILI, 1))

ILI_1112 <- ili_current %>%
  filter(season == "2011/2012") %>%
  select(season, location, week, ILI) %>%
  mutate(ILI = round(ILI, 1))

ILI_1213 <- ili_current %>%
  filter(season == "2012/2013") %>%
  select(season, location, week, ILI) %>%
  mutate(ILI = round(ILI, 1))

ILI_1314 <- readRDS("Data/ili_init_pub_list.rds")[['201428']] %>%
  filter(!location %in% state.name) %>%
  mutate(order_week = week_inorder(week, "2013/2014")) %>%
  arrange(location, order_week) %>%
  select(location, week, ILI) %>%
  mutate(ILI = round(ILI, 1))

ILI_1415 <- readRDS("Data/ili_init_pub_list.rds")[['201528']] %>%
  filter(!location %in% state.name) %>%
  mutate(order_week = week_inorder(week, "2014/2015")) %>%
  arrange(location, order_week) %>%
  select(location, week, ILI) %>%
  mutate(ILI = round(ILI, 1))

ILI_1516 <- readRDS("Data/ili_init_pub_list.rds")[['201628']] %>%
  filter(!location %in% state.name) %>%
  mutate(order_week = week_inorder(week, "2015/2016")) %>%
  arrange(location, order_week) %>%
  select(location, week, ILI) %>%
  mutate(ILI = round(ILI, 1))

ILI_1617 <- readRDS("Data/ili_init_pub_list.rds")[['201728']] %>%
  filter(!location %in% state.name) %>%
  mutate(order_week = week_inorder(week, "2016/2017")) %>%
  arrange(location, order_week) %>%
  select(location, week, ILI) %>%
  mutate(ILI = round(ILI, 1))

ILI_1718 <- readRDS("Data/ili_init_pub_list.rds")[['201828']] %>%
  filter(!location %in% state.name) %>%
  mutate(order_week = week_inorder(week, "2017/2018")) %>%
  arrange(location, order_week) %>%
  select(location, week, ILI) %>%
  mutate(ILI = round(ILI, 1))

ILI_1819 <- readRDS("Data/ili_init_pub_list.rds")[['201928']] %>%
  filter(!location %in% state.name) %>%
  mutate(order_week = week_inorder(week, "2018/2019")) %>%
  arrange(location, order_week) %>%
  select(location, week, ILI) %>%
  mutate(ILI = round(ILI, 1))

# Create truth ----------------------------------------------------------------
truth_1011 <- create_truth(fluview = FALSE, year = 2010, weekILI = ILI_1011,
                           challenge = "ilinet", start_wk = 40, end_wk = 20)

truth_1112 <- create_truth(fluview = FALSE, year = 2011, weekILI = ILI_1112,
                           challenge = "ilinet", start_wk = 40, end_wk = 20)

truth_1213 <- create_truth(fluview = FALSE, year = 2012, weekILI = ILI_1213,
                           challenge = "ilinet", start_wk = 40, end_wk = 20)

truth_1314 <- create_truth(fluview = FALSE, year = 2013, weekILI = ILI_1314,
                           challenge = "ilinet", start_wk = 40, end_wk = 20)

truth_1415 <- create_truth(fluview = FALSE, year = 2014, weekILI = ILI_1415,
                           challenge = "ilinet", start_wk = 40, end_wk = 20)

truth_1516 <- create_truth(fluview = FALSE, year = 2015, weekILI = ILI_1516,
                           challenge = "ilinet", start_wk = 40, end_wk = 20)

truth_1617 <- create_truth(fluview = FALSE, year = 2016, weekILI = ILI_1617,
                           challenge = "ilinet", start_wk = 40, end_wk = 20)

truth_1718 <- create_truth(fluview = FALSE, year = 2017, weekILI = ILI_1718,
                           challenge = "ilinet", start_wk = 40, end_wk = 20)

truth_1819 <- create_truth(fluview = FALSE, year = 2018, weekILI = ILI_1819,
                           challenge = "ilinet", start_wk = 40, end_wk = 20)


# Expand observed truth to include all bins that will be counted as correct ----
# exp_truth_1011 <- expand_truth(truth_1011, week_expand = 1, percent_expand = 5,
#                                challenge = "ilinet")
# 
# exp_truth_1112 <- expand_truth(truth_1112, week_expand = 1, percent_expand = 5,
#                                challenge = "ilinet")
# 
# exp_truth_1213 <- expand_truth(truth_1213, week_expand = 1, percent_expand = 5,
#                                challenge = "ilinet")
# 
# exp_truth_1314 <- expand_truth(truth_1314, week_expand = 1, percent_expand = 5,
#                                challenge = "ilinet")
# 
# exp_truth_1415 <- expand_truth(truth_1415, week_expand = 1, percent_expand = 5,
#                                challenge = "ilinet")
# 
# exp_truth_1516 <- expand_truth(truth_1516, week_expand = 1, percent_expand = 5,
#                                challenge = "ilinet")
# 
# exp_truth_1617 <- expand_truth(truth_1617, week_expand = 1, percent_expand = 5,
#                                challenge = "ilinet")
# 
# exp_truth_1718 <- expand_truth(truth_1718, week_expand = 1, percent_expand = 5,
#                                challenge = "ilinet")
# 
# exp_truth_1819 <- expand_truth(truth_1819, week_expand = 1, percent_expand = 5,
#                                challenge = "ilinet")

# Create evaluation period -----
eval_period_1011 <- create_eval_period(ILI_1011, truth_1011, "2010/2011")

eval_period_1112 <- create_eval_period(ILI_1112, truth_1112, "2011/2012")

eval_period_1213 <- create_eval_period(ILI_1213, truth_1213, "2012/2013")

eval_period_1314 <- create_eval_period(ILI_1314, truth_1314, "2013/2014")

eval_period_1415 <- create_eval_period(ILI_1415, truth_1415, "2014/2015")

eval_period_1516 <- create_eval_period(ILI_1516, truth_1516, "2015/2016")

eval_period_1617 <- create_eval_period(ILI_1617, truth_1617, "2016/2017")

eval_period_1718 <- create_eval_period(ILI_1718, truth_1718, "2017/2018")

eval_period_1819 <- create_eval_period(ILI_1819, truth_1819, "2018/2019")


# Score entries -----
full_scores_1011 <- calc_scores(forecasts_1011, truth_1011, 
                                season = "2010/2011", exclude = FALSE, 
                                eval = FALSE)

full_scores_1112 <- calc_scores(forecasts_1112, truth_1112, 
                                season = "2011/2012", exclude = FALSE, 
                                eval = FALSE)

full_scores_1213 <- calc_scores(forecasts_1213, truth_1213, 
                                season = "2012/2013", exclude = FALSE, 
                                eval = FALSE)

full_scores_1314 <- calc_scores(forecasts_1314, truth_1314, 
                                season = "2013/2014", exclude = FALSE, 
                                eval = FALSE)

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

eval_scores_1011 <- calc_scores(forecasts_1011, truth_1011, 
                                season = "2010/2011", exclude = FALSE, 
                                eval = TRUE, eval_period = eval_period_1011)

eval_scores_1112 <- calc_scores(forecasts_1112, truth_1112, 
                                season = "2011/2012", exclude = FALSE, 
                                eval = TRUE, eval_period = eval_period_1112)

eval_scores_1213 <- calc_scores(forecasts_1213, truth_1213, 
                                season = "2012/2013", exclude = FALSE, 
                                eval = TRUE, eval_period = eval_period_1213)

eval_scores_1314 <- calc_scores(forecasts_1314, truth_1314, 
                                season = "2013/2014", exclude = FALSE, 
                                eval = TRUE, eval_period = eval_period_1314)

eval_scores_1415 <- calc_scores(forecasts_1415, truth_1415, 
                                season = "2014/2015", exclude = FALSE, 
                                eval = TRUE, eval_period = eval_period_1415)

eval_scores_1516 <- calc_scores(forecasts_1516, truth_1516, 
                                season = "2015/2016", exclude = FALSE, 
                                eval = TRUE, eval_period = eval_period_1516)

eval_scores_1617 <- calc_scores(forecasts_1617, truth_1617, 
                                season = "2016/2017", exclude = FALSE, 
                                eval = TRUE, eval_period = eval_period_1617)

eval_scores_1718 <- calc_scores(forecasts_1718, truth_1718, 
                                season = "2017/2018", exclude = FALSE, 
                                eval = TRUE, eval_period = eval_period_1718)

eval_scores_1819 <- calc_scores(forecasts_1819, truth_1819, 
                                season = "2018/2019", exclude = FALSE, 
                                eval = TRUE, eval_period = eval_period_1819)

all_full_scores <- bind_rows(full_scores_1011, full_scores_1112, 
                             full_scores_1213, full_scores_1314,
                             full_scores_1415, full_scores_1516, 
                             full_scores_1617, full_scores_1718,
                             full_scores_1819)

all_eval_scores <- bind_rows(eval_scores_1011, eval_scores_1112, 
                             eval_scores_1213, eval_scores_1314,
                             eval_scores_1415, eval_scores_1516, 
                             eval_scores_1617, eval_scores_1718,
                             eval_scores_1819)

# Save scores -----
saveRDS(all_full_scores, "Data/full_model_scores.Rds")
saveRDS(all_eval_scores, "Data/eval_model_scores.Rds")

