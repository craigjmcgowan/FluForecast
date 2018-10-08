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
#   inner_join(select(virologic_combined, location, season, year, week, h1_per_samples,
#                    h3_per_samples, b_per_samples),
#             by = c("location", "season", "year", "week")) %>%
#   inner_join(gtrend_US_flu_merge,
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
#   # Add Google Trends data by location
#   left_join(bind_rows(mutate(gtrend_US_flu_merge, location = "US National"),
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
#   # Replace missing backfill data with random draw from same season/week observed values
#   mutate(
#     sim_backfill = map2_dbl(location, week,
#                         ~ sample(ili_backfill$backfill[ili_backfill$location == .x &
#                                                          ili_backfill$week == .y], 1)),
#     backfill = case_when(
#       !is.na(backfill) ~ backfill,
#       TRUE ~ sim_backfill
#       )
#     ) %>%
#   select(-sim_backfill)
# 
# save(nested_truth, flu_data_merge, file = "Data/truth_and_data.Rdata")

### Decisions to make ###

## Transformation of data
## Number of Fourier terms
## Structure of ARIMA errors
## What virologic/Gtrend data to include
## If/how to include estimates of backfill\

##### Transform data #####
load("Data/transform_fits.Rdata")
# transform_model_fit_data <- tibble(season = c("2010/2011", "2011/2012", "2012/2013", "2013/2014",
#                                 "2015/2016", "2016/2017", "2017/2018")) %>%
#   mutate(train_data = map(season,
#                           ~ filter(flu_data_merge, year <= as.numeric(substr(., 6, 9)),
#                                    season != paste0(substr(., 6, 9), "/",
#                                                     as.numeric(substr(., 6, 9)) + 1),
#                                    season != .))) %>%
#   unnest() %>%
#   # Nest by season and location
#   nest(-season, -location) %>%
#   # Create time series of ILI
#   mutate(data = map(data,
#                    ~ mutate(.x,
#                             ILI = ts(ILI, frequency = 52, start = c(2006, 40)))),
#          group = rep(1:length(cluster), length.out = nrow(.)))
# 
# transform_fits_by_group <- transform_model_fit_data %>%
#   partition(group, cluster = cluster)
# 
# transform_fits_by_group %>% 
#   cluster_library(c("tidyverse", "forecast", "lubridate", "FluSight", "MMWRweek"))
# 
# transform_model_fits <- transform_fits_by_group %>%
#   mutate(no_trans = map(data,
#                      ~ auto.arima(.$ILI, xreg = fourier(.$ILI, K = 1),
#                                   seasonal = FALSE)),
#          log = map(data,
#                      ~ auto.arima(.$ILI, xreg = fourier(.$ILI, K = 1),
#                                   seasonal = FALSE, lambda = 0)),
#          Box_Cox = map(data,
#                      ~ auto.arima(.$ILI, xreg = fourier(.$ILI, K = 1),
#                                   seasonal = FALSE, lambda = "auto"))) %>%
#   collect() %>%
#   as.tibble() %>%
#   ungroup() %>%
#   select(-data) %>%
#   gather(key = "model", value = "fit", no_trans, log, Box_Cox)
# 
# save(transform_model_fits, file = "Data/transform_fits.Rdata")

# Set up data for model fitting in parallel
transform_model_data <- crossing(model = c("no_trans", "log", "Box_Cox"),
                               season = c("2010/2011", "2011/2012", "2012/2013", "2013/2014",
                                          "2014/2015", "2015/2016", "2016/2017", "2017/2018"),
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
  left_join(transform_model_fits, by = c("season", "location", "model")) %>%
  # Set up grouping for parallel
  mutate(group = rep(1:length(cluster), length.out = nrow(.)))


# Set up party_df and load necessary libraries and functions 
transform_by_group <- transform_model_data %>%
  partition(group, cluster = cluster)

transform_by_group %>% 
  cluster_library(c("tidyverse", "forecast", "lubridate", "FluSight", "MMWRweek")) %>% 
  cluster_assign_value("flu_data_merge", flu_data_merge) %>%
  cluster_assign_value("ili_init_pub_list", ili_init_pub_list) %>%
  cluster_assign_value("fit_to_forecast", fit_to_forecast) %>%
  cluster_assign_value("sample_predictive_trajectories_arima", 
                       sample_predictive_trajectories_arima)

# Create forecasts for different transformations
transform_forecasts <- transform_by_group %>%
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
    xreg = map2(pred_data, model,
                ~ fourier(.x$ILI, K = 1)),
    max_week = ifelse(season == "2014/2015", 53, 52),
    # Fit models
    pred_fit = pmap(list(pred_data, fit, xreg),
                    ~ Arima(..1$ILI, xreg = ..3, model = ..2)),
    pred_results = pmap(
      list(pred_fit, pred_data, location, season, model, max_week),
      ~ fit_to_forecast(object = ..1,
                        xreg = fourier(..2$ILI, K = 1,
                                       h = ..6 - 17 - nrow(..2[..2$season == ..4, ])),
                        pred_data = ..2,
                        location = ..3,
                        season = ..4,
                        max_week = ..6,
                        npaths = 250)
    )
  ) %>%
  collect() %>%
  as.tibble()

transform_forecasts <- transform_forecasts %>% ungroup() %>%
  select(season, model, location, week, pred_results)

save(transform_forecasts, file = "Data/transform_forecasts.Rdata") 

# Normalize probabilities and score forecasts 
transform_scores <- transform_forecasts %>% 
  mutate(pred_results = map2(pred_results, location,
                             ~ mutate(.x, location = .y) %>%
                               normalize_probs())) %>%
  select(season, model, week, pred_results) %>%
  unnest() %>%
  nest(-season, -model, -week) %>%
  left_join(nested_truth, by = "season") %>%
  mutate(scores = map2(data, exp_truth,
                       ~ score_entry(.x, .y)),
         eval_scores = pmap(list(scores, eval_period, season),
                            ~ create_eval_scores(..1, ..2, ..3))) %>%
  select(season, model, eval_scores) %>%
  unnest() 

save(transform_scores, file = "Data/CV_Transform_Scores.Rdata")

# Determine best K value for each region
best_transform_cv <- transform_scores %>%
  group_by(location, model) %>%
  summarize(avg_score = mean(score)) %>%
  filter(avg_score == max(avg_score)) %>%
  ungroup() %>%
  select(location, transform = model)

save(best_transform_cv, file = "Data/CV_Transform_terms.Rdata")

##### Number of Fourier terms #####
load("Data/fourier_fits.Rdata")
load("Data/CV_Transform_terms.Rdata")
# fourier_model_fits <- tibble(season = c("2010/2011", "2011/2012", "2012/2013", "2013/2014",
#                                 "2015/2016", "2016/2017", "2017/2018")) %>%
#   mutate(train_data = map(season,
#                           ~ filter(flu_data_merge, year <= as.numeric(substr(., 6, 9)),
#                                    season != paste0(substr(., 6, 9), "/",
#                                                     as.numeric(substr(., 6, 9)) + 1),
#                                    season != .))) %>%
#   unnest() %>%
#   # Nest by season and location
#   nest(-season, -location) %>%
#   # Join lambda value for arima model
#   left_join(best_transform_cv, by = "location") %>%
#   # Create time series of ILI and numerical value of lambda
#   mutate(data = map(data,
#                    ~ mutate(.x,
#                             ILI = ts(ILI, frequency = 52, start = c(2006, 40)))),
#          lambda = unlist(map(data,
#                              ~ BoxCox.lambda(.$ILI))),
#          lambda = case_when(transform == "no_trans" ~ NA_real_,
#                             transform == "log" ~ 0,
#                             transform == "Box_Cox" ~ lambda)) %>%
#   # Fit model
#   mutate(fit_1 = map2(data, lambda,
#                      ~ auto.arima(.x$ILI, xreg = fourier(.x$ILI, K = 1),
#                                   seasonal = FALSE, lambda = .y)),
#          fit_2 = map2(data, lambda,
#                       ~ auto.arima(.x$ILI, xreg = fourier(.x$ILI, K = 2),
#                                    seasonal = FALSE, lambda = .y)),
#          fit_3 = map2(data, lambda,
#                       ~ auto.arima(.x$ILI, xreg = fourier(.x$ILI, K = 3),
#                                    seasonal = FALSE, lambda = .y)),
#          fit_4 = map2(data, lambda,
#                       ~ auto.arima(.x$ILI, xreg = fourier(.x$ILI, K = 4),
#                                    seasonal = FALSE, lambda = .y)),
#          fit_5 = map2(data, lambda,
#                       ~ auto.arima(.x$ILI, xreg = fourier(.x$ILI, K = 5),
#                                    seasonal = FALSE, lambda = .y)),
#          fit_6 = map2(data, lambda,
#                       ~ auto.arima(.x$ILI, xreg = fourier(.x$ILI, K = 6),
#                                    seasonal = FALSE, lambda = .y)),
#          fit_7 = map2(data, lambda,
#                       ~ auto.arima(.x$ILI, xreg = fourier(.x$ILI, K = 7),
#                                    seasonal = FALSE, lambda = .y)),
#          fit_8 = map2(data, lambda,
#                       ~ auto.arima(.x$ILI, xreg = fourier(.x$ILI, K = 8),
#                                    seasonal = FALSE, lambda = .y)),
#          fit_9 = map2(data, lambda,
#                       ~ auto.arima(.x$ILI, xreg = fourier(.x$ILI, K = 9),
#                                    seasonal = FALSE, lambda = .y)),
#          fit_10 = map2(data, lambda,
#                        ~ auto.arima(.x$ILI, xreg = fourier(.x$ILI, K = 10),
#                                     seasonal = FALSE, lambda = .y)),
#          fit_11 = map2(data, lambda,
#                        ~ auto.arima(.x$ILI, xreg = fourier(.x$ILI, K = 11),
#                                     seasonal = FALSE, lambda = .y)),
#          fit_12 = map2(data, lambda,
#                        ~ auto.arima(.x$ILI, xreg = fourier(.x$ILI, K = 12),
#                                     seasonal = FALSE, lambda = .y))) %>%
#   select(-data) %>%
#   gather(key = "model", value = "fit", fit_1:fit_12)
# 
# save(fourier_model_fits, file = "Data/fourier_fits.Rdata")

# Set up data for model fitting in parallel
fourier_model_data <- crossing(model = c("fit_1", "fit_2", "fit_3", "fit_4",
                                         "fit_5", "fit_6", "fit_7", "fit_8", 
                                         "fit_9", "fit_10", "fit_11", "fit_12"),
                       season = c("2010/2011", "2011/2012", "2012/2013", "2013/2014",
                                  "2014/2015", "2015/2016", "2016/2017", "2017/2018"),
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
  left_join(fourier_model_fits, by = c("season", "location", "model")) %>%
  # Set up grouping for parallel
  mutate(group = rep(1:4, length.out = nrow(.)))


# Set up party_df and load necessary libraries and functions
# load("Data/fourier_scores.Rdata")
fourier_scores <- tibble()
start_time <- Sys.time()
for (this_season in c("2010/2011", "2011/2012", "2012/2013", "2013/2014", 
                      "2014/2015", "2015/2016", "2016/2017", "2017/2018")) {

  fourier_cluster <- create_cluster(cores = parallel::detectCores())
  
  fourier_by_group <- fourier_model_data %>%
    filter(season == this_season) %>%
    partition(group, cluster = fourier_cluster)
  
  fourier_by_group %>% 
    cluster_library(c("tidyverse", "forecast", "lubridate", "FluSight", "MMWRweek")) %>% 
    cluster_assign_value("flu_data_merge", flu_data_merge) %>%
    cluster_assign_value("ili_init_pub_list", ili_init_pub_list) %>%
    cluster_assign_value("fit_to_forecast", fit_to_forecast) %>%
    cluster_assign_value("sample_predictive_trajectories_arima", 
                         sample_predictive_trajectories_arima)

  # Create forecasts for Fourier terms
  fourier_forecasts_temp <- fourier_by_group %>%
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
      # Set up Fourier data for forecasting
      xreg = map2(pred_data, model,
                  ~ fourier(.x$ILI, K = as.numeric(str_extract(.y, "[1-9][0-9]|[0-9]")))),
      max_week = ifelse(season == "2014/2015", 53, 52),
      # Fit models
      pred_fit = pmap(list(pred_data, fit, xreg),
                      ~ Arima(..1$ILI, xreg = ..3, model = ..2)),
      pred_results = pmap(
        list(pred_fit, pred_data, location, season, model, max_week),
        ~ fit_to_forecast(object = ..1,
                          xreg = fourier(..2$ILI, K = as.numeric(str_extract(..5, "[1-9][0-9]|[0-9]")),
                                         h = ..6 - 17 - nrow(..2[..2$season == ..4, ])),
                          pred_data = ..2,
                          location = ..3,
                          season = ..4,
                          max_week = ..6,
                          npaths = 250)
        )
      ) %>%
    collect() %>%
    as.tibble()
    
  fourier_forecasts_temp <- fourier_forecasts_temp %>% ungroup() %>%
    select(season, model, location, week, pred_results) 
  
  save(fourier_forecasts_temp, 
       file = paste0("Data/fourier_forecasts_", 
                                             substr(this_season, 1, 4), ".Rdata")) 

  # Normalize probabilities and score forecasts 
  fourier_scores <- fourier_forecasts_temp %>% 
    mutate(pred_results = map2(pred_results, location,
                               ~ mutate(.x, location = .y) %>%
                                 normalize_probs())) %>%
    select(season, model, week, pred_results) %>%
    unnest() %>%
    nest(-season, -model, -week) %>%
    left_join(nested_truth, by = "season") %>%
    mutate(scores = map2(data, exp_truth,
                         ~ score_entry(.x, .y)),
           eval_scores = pmap(list(scores, eval_period, season),
                              ~ create_eval_scores(..1, ..2, ..3))) %>%
    select(season, model, eval_scores) %>%
    unnest() %>%
    bind_rows(fourier_scores)
  
  save(fourier_scores, file = "Data/CV_Fourier_Scores.Rdata")
  
  rm(fourier_forecasts_temp, fourier_by_group, fourier_cluster)
}

# Determine best K value for each region
best_k_cv <- fourier_scores %>%
  group_by(location, model) %>%
  summarize(avg_score = mean(score)) %>%
  filter(avg_score == max(avg_score)) %>%
  mutate(K = as.numeric(str_extract(model, "[1-9][0-9]|[0-9]"))) %>%
  ungroup() %>%
  select(location, K)

save(best_k_cv, file = "Data/CV_Fourier_terms.Rdata")

##### ARIMA structure for error terms #####

load("Data/CV_Fourier_terms.Rdata")
load("Data/CV_Transform_terms.Rdata")
load("Data/arima_fits.Rdata")
# arima_model_fit_data <- crossing(season = c("2010/2011", "2011/2012", "2012/2013",
#                                             "2013/2014", "2014/2015", "2015/2016",
#                                             "2016/2017", "2017/2018"),
#                                 arima_1 = 0:3,
#                                 arima_2 = 0:1,
#                                 arima_3 = 0:3) %>%
#   mutate(train_data = map(season,
#                           ~ filter(flu_data_merge, year <= as.numeric(substr(., 6, 9)),
#                                    season != paste0(substr(., 6, 9), "/",
#                                                     as.numeric(substr(., 6, 9)) + 1),
#                                    season != .))) %>%
#   unnest() %>%
#   # Nest by season and location
#   nest(-season, -location, -arima_1, -arima_2, -arima_3) %>%
#   # Create time series of ILI
#   mutate(data = map(data,
#                    ~ mutate(.x,
#                             ILI = ts(ILI, frequency = 52, start = c(2006, 40))))) %>%
#   # Join lambda and fourier values
#   left_join(best_transform_cv, by = "location") %>%
#   left_join(best_k_cv, by = "location") %>%
#   # Set up for parallel
#   mutate(group = rep(1:length(cluster), length.out = nrow(.)),
#          lambda = unlist(map(data,
#                              ~ BoxCox.lambda(.$ILI))),
#          lambda = case_when(model == "no_trans" ~ NA_real_,
#                             model == "log" ~ 0,
#                             model == "Box_Cox" ~ lambda)) %>%
#   select(-model)
# 
# # Set up clusters
# arima_fit_parallel <- arima_model_fit_data %>%
#   partition(group, cluster = cluster)
# arima_fit_parallel %>%
#   cluster_library(c("tidyverse", "forecast", "lubridate", "FluSight", "MMWRweek"))
# 
# # Fit ARIMA models
# arima_model_fits <- arima_fit_parallel %>%
#   mutate(fit = pmap(
#     list(data, arima_1, arima_2, arima_3, K, lambda),
#     ~ try(Arima(..1$ILI, order = c(..2, ..3, ..4),
#             xreg = fourier(..1$ILI, K = ..5),
#             lambda = ..6), silent = TRUE)
#     )) %>%
#   select(-data) %>%
#   collect() %>%
#   as.tibble() %>%
#   ungroup() %>%
#   filter(!grepl("Error", fit))
# 
# arima_model_fits <- arima_model_fits%>%
#   filter(!grepl("Error", fit))
# 
# save(arima_model_fits, file = "Data/arima_fits.Rdata")

# Set up data for forecast creation in parallel
arima_model_data_setup <- crossing(season = c("2010/2011", "2011/2012", "2012/2013",
                                              "2013/2014", "2014/2015", "2015/2016",
                                              "2016/2017", "2017/2018"),
                             arima_1 = 0:3,
                             arima_2 = 0:1,
                             arima_3 = 0:3,
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
  right_join(select(arima_model_fits, -group), 
            by = c("location", "season", "arima_1", "arima_2", "arima_3")) %>%
  # Set up grouping for parallel
  mutate(group = rep(1:4, length.out = nrow(.)))
  

# Create and score forecasts in single-season parallel chunks
# Clear forecasts from memory after each season

start_time <- Sys.time()
for(this_season in unique(arima_model_data_setup$season)) {
  
  for_cluster <- create_cluster(cores = parallel::detectCores())

  arima_model_data_parallel <- arima_model_data_setup %>%
    filter(season == this_season) %>%
    partition(group, cluster = for_cluster)
  
  arima_model_data_parallel %>%
    cluster_library(c("tidyverse", "forecast", "lubridate", "FluSight", "MMWRweek")) %>% 
    cluster_assign_value("flu_data_merge", flu_data_merge) %>%
    cluster_assign_value("ili_init_pub_list", ili_init_pub_list) %>%
    cluster_assign_value("fit_to_forecast", fit_to_forecast) %>%
    cluster_assign_value("sample_predictive_trajectories_arima", 
                         sample_predictive_trajectories_arima)
  

  # Create forecasts for different ARIMA structures
  arima_forecasts <- arima_model_data_parallel %>%
    mutate(pred_data = pmap(list(season, week, location, epiweek), 
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
                       mutate(ILI = ts(ILI, frequency = 52, start = c(2006, 40))))) %>%
    # Set up Fourier data for forecasting
    mutate(xreg = map2(pred_data, K,
                       ~ fourier(.x$ILI, K = .y)),
           max_week = ifelse(season == "2014/2015", 53, 52)) %>%
    # Fit ARIMA model based on most recent data using prior fit
    mutate(pred_fit = pmap(list(pred_data, fit, xreg),
                           ~ Arima(..1$ILI, xreg = ..3, model = ..2))) %>%
    # Create predicted results
    mutate(pred_results = pmap(
      list(pred_fit, pred_data, location, season, max_week, K),
      ~ fit_to_forecast(object = ..1,
                        xreg = fourier(..2$ILI, K = ..6,
                                       h = ..5 - 17 - nrow(..2[..2$season == ..4, ])),
                        pred_data = ..2,
                        location = ..3,
                        season = ..4,
                        max_week = ..5,
                        npaths = 250)
      )) %>%
    collect() %>%
    as.tibble() %>%
    ungroup() %>%
    select(season, arima_1:arima_3, location, week, pred_results)
 
  # Normalize probabilities and score forecasts 
  arima_scores <- arima_forecasts %>%
    mutate(pred_results = map2(pred_results, location,
                               ~ mutate(.x, location = .y) %>%
                                 normalize_probs())) %>%
    select(season, arima_1:arima_3, week, pred_results) %>%
    unnest() %>%
    nest(-season, -arima_1, -arima_2, -arima_3, -week) %>%
    left_join(nested_truth, by = "season") %>%
    mutate(scores = map2(data, exp_truth,
                         ~ score_entry(.x, .y)),
           eval_scores = pmap(list(scores, eval_period, season),
                              ~ create_eval_scores(..1, ..2, ..3))) %>%
    select(season, arima_1:arima_3, eval_scores) %>%
    unnest() 
  
  save(arima_scores, file = paste0("Data/CV_ARIMA_Scores_HHS2_", substr(this_season, 1, 4),
                                   ".Rdata"))
  
  rm(for_cluster, arima_forecasts, arima_model_data_parallel)

}
Sys.time() - start_time

load("Data/CV_ARIMA_Scores_2010.Rdata")
arima_scores_1011 <- arima_scores
load("Data/CV_ARIMA_Scores_2011.Rdata")
arima_scores_1112 <- arima_scores
load("Data/CV_ARIMA_Scores_2012.Rdata")
arima_scores_1213 <- arima_scores
load("Data/CV_ARIMA_Scores_2013.Rdata")
arima_scores_1314 <- arima_scores
load("Data/CV_ARIMA_Scores_2014.Rdata")
arima_scores_1415 <- arima_scores
load("Data/CV_ARIMA_Scores_2015.Rdata")
arima_scores_1516 <- arima_scores
load("Data/CV_ARIMA_Scores_2016.Rdata")
arima_scores_1617 <- arima_scores
load("Data/CV_ARIMA_Scores_2017.Rdata")
arima_scores_1718 <- arima_scores

# Determine best ARIMA model for each region
best_arima_cv <- bind_rows(arima_scores_1011, arima_scores_1112,
                           arima_scores_1213, arima_scores_1314,
                           arima_scores_1415, arima_scores_1516,
                           arima_scores_1617, arima_scores_1718)%>%
  group_by(location, arima_1, arima_2, arima_3) %>%
  summarize(avg_score = mean(score)) %>%
  group_by(location) %>%
  filter(avg_score == max(avg_score)) %>%
  ungroup() %>%
  select(location, arima_1:arima_3, avg_score)

save(best_arima_cv, file = "Data/CV_ARIMA_terms.Rdata")

##### Additional covariates #####
load("Data/CV_Transform_terms.Rdata")
load("data/CV_Fourier_terms.Rdata")
load("Data/CV_ARIMA_terms.Rdata")

# Models to test:
# No covariates
# Backfill only
# National Google Trends data only
# Regional GTrends data only
# Flu % only
# National Google Trends data and backfill
# Regional Google Trends data and backfill
# Flu % and backfill
# National Gtrends and flu %
# Regional Gtrends and flu %
# National Gtrends, flu, and backfill
# Regional Gtrends, flu, and backfill

load("Data/covar_fits.Rdata")
# covar_model_fits <- crossing(season = c("2010/2011", "2011/2012", "2012/2013",
#                                         "2013/2014", "2014/2015", "2015/2016",
#                                         "2016/2017", "2017/2018"),
#                              model = c("ARIMA only", "National Gtrends", "FluVirus", 
#                                        "Backfill", "National Gtrends & Backfill",
#                                        "FluVirus & Backfill",
#                                        "National Gtrends & FluVirus",
#                                        "National Gtrends, FluVirus, & Backfill",
#                                        "Regional Gtrends", "Regional Gtrends & Backfill",
#                                        "Regional Gtrends & FluVirus",
#                                        "Regional Gtrends, FluVirus, & Backfill")) %>%
#   mutate(train_data = map(season,
#                           ~ filter(flu_data_merge, year <= as.numeric(substr(., 6, 9)),
#                                    season != paste0(substr(., 6, 9), "/",
#                                                     as.numeric(substr(., 6, 9)) + 1),
#                                    season != .))) %>%
#   unnest() %>%
#   # Nest by season and location
#   nest(-season, -location, -model) %>%
#   # Create time series of ILI
#   mutate(data = map(data,
#                    ~ mutate(.x,
#                             ILI = ts(ILI, frequency = 52, start = c(2006, 40))))) %>%
#   # Merge best lambda, Fourier K value, and ARIMA structure by location
#   left_join(rename(best_transform_cv, transform = model), by = "location") %>%
#   mutate(lambda = unlist(map(data,
#                              ~ BoxCox.lambda(.$ILI))),
#          lambda = case_when(transform == "no_trans" ~ NA_real_,
#                             transform == "log" ~ 0,
#                             transform == "Box_Cox" ~ lambda)) %>%
#   left_join(best_k_cv, by = "location") %>%
#   left_join(best_arima_cv, by = "location") %>%
#   # Fit models
#   mutate(fit = case_when(
#     model == "ARIMA only" ~
#       pmap(
#         list(data, arima_1, arima_2, arima_3, K, lambda),
#         ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
#                 xreg = fourier(..1$ILI, K = ..5),
#                 lambda = ..6)
#       ),
#     model == "National Gtrends" ~
#       pmap(
#         list(data, arima_1, arima_2, arima_3, K, lambda),
#         ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
#                 xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
#                              ..1$hits),
#                 lambda = ..6)
#       ),
#     model == "Regional Gtrends" ~
#       pmap(
#         list(data, arima_1, arima_2, arima_3, K, lambda),
#         ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
#                 xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
#                              ..1$region_hits),
#                 lambda = ..6)
#       ),
#     model == "FluVirus" ~
#       pmap(
#         list(data, arima_1, arima_2, arima_3, K, lambda),
#         ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
#                 xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
#                              ..1$h1_per_samples, ..1$h3_per_samples,
#                              ..1$b_per_samples),
#                 lambda = ..6)
#       ),
#     model == "Backfill" ~ 
#       pmap(
#         list(data, arima_1, arima_2, arima_3, K, lambda),
#         ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
#                 xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
#                              ..1$backfill),
#                 lambda = ..6)
#       ),
#     model == "National Gtrends & Backfill" ~
#       pmap(
#         list(data, arima_1, arima_2, arima_3, K, lambda),
#         ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
#                 xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
#                              ..1$hits, ..1$backfill),
#                 lambda = ..6)
#       ),
#     model == "Regional Gtrends & Backfill" ~
#       pmap(
#         list(data, arima_1, arima_2, arima_3, K, lambda),
#         ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
#                 xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
#                              ..1$region_hits, ..1$backfill),
#                 lambda = ..6)
#       ),
#     model == "FluVirus & Backfill" ~
#       pmap(
#         list(data, arima_1, arima_2, arima_3, K, lambda),
#         ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
#                 xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
#                              ..1$h1_per_samples, ..1$h3_per_samples,
#                              ..1$b_per_samples, ..1$backfill),
#                 lambda = ..6)
#       ),
#     model == "National Gtrends & FluVirus" ~
#       pmap(
#         list(data, arima_1, arima_2, arima_3, K, lambda),
#         ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
#                 xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
#                              ..1$hits, ..1$h1_per_samples,
#                              ..1$h3_per_samples, ..1$b_per_samples),
#                 lambda = ..6)
#       ),
#     model == "Regional Gtrends & FluVirus" ~
#       pmap(
#         list(data, arima_1, arima_2, arima_3, K, lambda),
#         ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
#                 xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
#                              ..1$region_hits, ..1$h1_per_samples,
#                              ..1$h3_per_samples, ..1$b_per_samples),
#                 lambda = ..6)
#       ),
#     model == "National Gtrends, FluVirus, & Backfill" ~
#       pmap(
#         list(data, arima_1, arima_2, arima_3, K, lambda),
#         ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
#                 xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
#                              ..1$hits, ..1$h1_per_samples,
#                              ..1$h3_per_samples, ..1$b_per_samples,
#                              ..1$backfill),
#                 lambda = ..6)
#       ),
#     model == "Regional Gtrends, FluVirus, & Backfill" ~
#       pmap(
#         list(data, arima_1, arima_2, arima_3, K, lambda),
#         ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
#                 xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
#                              ..1$region_hits, ..1$h1_per_samples,
#                              ..1$h3_per_samples, ..1$b_per_samples,
#                              ..1$backfill),
#                 lambda = ..6)
#       )
#     ))
# 
# save(covar_model_fits, file = "Data/covar_fits.Rdata")

# Set up data for forecast creation in parallel
covar_model_data_setup <- crossing(season = c("2010/2011", "2011/2012", "2012/2013",
                                              "2013/2014", "2014/2015", "2015/2016",
                                              "2016/2017", "2017/2018"),
                                   model = c("ARIMA only", "National Gtrends", "FluVirus", 
                                             "Backfill", "National Gtrends & Backfill",
                                             "FluVirus & Backfill",
                                             "National Gtrends & FluVirus",
                                             "National Gtrends, FluVirus, & Backfill",
                                             "Regional Gtrends", "Regional Gtrends & Backfill",
                                             "Regional Gtrends & FluVirus",
                                             "Regional Gtrends, FluVirus, & Backfill"),
                                   week = c(43:71),
                                   location = unique(flu_data_merge$location)) %>%
  filter(week < 71 | season == "2014/2015") %>%
  mutate(
    epiweek = case_when(
      season == "2014/2015" & week > 53 ~ 
        as.numeric(paste0(substr(season, 6, 9), 
                          str_pad(week - 53, 2, "left", "0"))),
      week > 52 ~ 
        as.numeric(paste0(substr(season, 6, 9), 
                          str_pad(week - 52, 2, "left", "0"))),
      TRUE ~ as.numeric(paste0(substr(season, 1, 4), str_pad(week, 2, "left", "0")))
    ),
    prev_season = case_when(
      season == "2010/2011" ~ "2007/2008",
      TRUE ~ paste0(as.numeric(substr(season, 1, 4)) - 1, "/",
                    substr(season, 1, 4))
    )
  ) %>%
  inner_join(covar_model_fits, 
            by = c("location", "season", "model")) %>%
  # Set up grouping for parallel
  mutate(group = rep(1:4, length.out = nrow(.)))

start_time <- Sys.time()
for (this_season in unique(covar_model_data_setup$season)) {
  
  cluster <- create_cluster(cores = parallel::detectCores())
  
  # Set up dataset for forecasting in parallel
  covar_model_parallel <- covar_model_data_setup %>%
    filter(season == this_season) %>%
    partition(group, cluster = cluster)
  
  covar_model_parallel %>%
    cluster_library(c("tidyverse", "forecast", "lubridate", "FluSight", "MMWRweek")) %>% 
    cluster_assign_value("flu_data_merge", flu_data_merge) %>%
    cluster_assign_value("ili_init_pub_list", ili_init_pub_list) %>%
    cluster_assign_value("fit_to_forecast", fit_to_forecast) %>%
    cluster_assign_value("sample_predictive_trajectories_arima", 
                         sample_predictive_trajectories_arima)
  
  # Create forecasts in parallel
  # For Gtrends - take next observed value since it would have been known at the time
  #   and for future values use seasonal naive adjusted by ratio of current year's 1 wk
  #   ahead value to the same value 1 year ago
  # For virus percentages, just carry forward most recent cumulative percentage into 
  #   the future. This will be unstable early in the season but will become more
  #   stable and informative as the season goes on
  # For backfill - use predicted value from linear model of previous backfits. For
  #   seasons before 2014/2015, use simulated backfills in 2004/2005-2007/2008 seasons
  #   to help inform predictions. For 2014/2015 onwards, only use observed backfills
  #   in order to eliminate impact of those simulated seasons
  
  covar_forecasts <- covar_model_parallel %>%
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
      # Set up Fourier data for forecasting
      xreg = case_when(
        model == "ARIMA only" ~ 
          map2(pred_data, K,
               ~ as.data.frame(fourier(.x$ILI, K = .y))),
        model == "National Gtrends" ~
          pmap(list(pred_data, K), 
               ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2)),
                       ..1$hits)),
        model == "Regional Gtrends" ~
          pmap(list(pred_data, K), 
               ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2)),
                       ..1$region_hits)),
        model == "FluVirus" ~ 
          pmap(list(pred_data, K),
               ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2)),
                       ..1$h1_per_samples, ..1$h3_per_samples,
                       ..1$b_per_samples)),
        model == "Backfill" ~ 
          pmap(list(pred_data, K),
               ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2)),
                       ..1$backfill)),
        model == "National Gtrends & Backfill" ~
          pmap(list(pred_data, K), 
               ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2)),
                       ..1$hits, ..1$backfill)),
        model == "Regional Gtrends & Backfill" ~
          pmap(list(pred_data, K), 
               ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2)),
                       ..1$region_hits, ..1$backfill)),
        model == "FluVirus & Backfill" ~ 
          pmap(list(pred_data, K),
               ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2)),
                       ..1$h1_per_samples, ..1$h3_per_samples,
                       ..1$b_per_samples, ..1$backfill)),
        model == "National Gtrends & FluVirus" ~ 
          pmap(list(pred_data, K),
               ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2)),
                     ..1$hits, ..1$h1_per_samples, 
                     ..1$h3_per_samples, ..1$b_per_samples)),
        model == "Regional Gtrends & FluVirus" ~ 
          pmap(list(pred_data, K),
               ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2)),
                       ..1$region_hits, ..1$h1_per_samples, 
                       ..1$h3_per_samples, ..1$b_per_samples)),
        model == "National Gtrends, FluVirus, & Backfill" ~ 
          pmap(list(pred_data, K),
               ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2)),
                       ..1$hits, ..1$h1_per_samples, 
                       ..1$h3_per_samples, ..1$b_per_samples,
                       ..1$backfill)),
        model == "Regional Gtrends, FluVirus, & Backfill" ~ 
          pmap(list(pred_data, K),
               ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2)),
                       ..1$region_hits, ..1$h1_per_samples, 
                       ..1$h3_per_samples, ..1$b_per_samples,
                       ..1$backfill))
        ),
      max_week = ifelse(season == "2014/2015", 53, 52),
      # Fit models
      pred_fit = pmap(list(pred_data, fit, xreg),
                      ~ Arima(..1$ILI, xreg = ..3, model = ..2)),
      # Create data frame of xreg terms for forecasting
      gtrend_forecast = pmap(
        list(pred_data, location, season, prev_season, week, max_week),
        ~ tibble(hits = c(flu_data_merge %>%
                            filter(location == ..2, season == ..3, 
                                   order_week == ..5 + 1) %>%
                            pull(hits),
                          flu_data_merge %>%
                            filter(location == ..2, season == ..4, 
                                   order_week > ..5 + 1, 
                                   order_week < ..6 + 26) %>%
                            pull(hits) *
                            mean(flu_data_merge %>%
                              filter(location == ..2, season == ..3,
                                     order_week <= ..5 + 1) %>%
                              pull(hits) /
                              flu_data_merge %>%
                                filter(location == ..2, season == ..4,
                                       order_week <= ..5 + 1) %>%
                                pull(hits))))
      ),
      reg_gtrend_forecast = pmap(
        list(pred_data, location, season, prev_season, week, max_week),
        ~ tibble(region_hits = c(flu_data_merge %>%
                            filter(location == ..2, season == ..3, 
                                   order_week == ..5 + 1) %>%
                            pull(region_hits),
                          flu_data_merge %>%
                            filter(location == ..2, season == ..4, 
                                   order_week > ..5 + 1, 
                                   order_week < ..6 + 26) %>%
                            pull(region_hits) *
                            mean(flu_data_merge %>%
                                   filter(location == ..2, season == ..3,
                                          order_week <= ..5 + 1) %>%
                                   pull(region_hits) /
                                   flu_data_merge %>%
                                   filter(location == ..2, season == ..4,
                                          order_week <= ..5 + 1) %>%
                                   pull(region_hits))))
      ),
      h1_per_forecast = pmap(
        list(pred_data, season, max_week),
        ~ rep(last(..1$h1_per_samples), ..3 - 14 - 
                nrow(..1[..1$season == ..2, ]))
      ),
      h3_per_forecast = pmap(
        list(pred_data, season, max_week),
        ~ rep(last(..1$h3_per_samples), ..3 - 14 - 
                nrow(..1[..1$season == ..2, ]))
      ),
      b_per_forecast = pmap(
        list(pred_data, season, max_week),
        ~ rep(last(..1$b_per_samples), ..3 - 14 - 
                nrow(..1[..1$season == ..2, ]))
      ),
      backfill = pmap(
        list(data, week, max_week, season),
        function(data, week, max_week, season) {
          out <- numeric()
          for(i in week:(max_week + 24)) {
            if(season %in% c("2010/2011", "2011/2012", "2012/2013")) {
              temp <- data[data$order_week == i, ]
            } else{
              temp <- data[data$order_week == i &
                             !data$season %in% c("2010/2011", "2011/2012", "2012/2013"), ]
            }
            
            out <- c(out, predict(lm(backfill ~ year, 
                                     data = data[data$order_week == i, ]),
                                  data.frame(year = as.numeric(substr(season, 1, 4)))))
          }
          out
        } 
      ),
      forecast_xreg = case_when(
        model == "ARIMA only" ~ 
          pmap(list(pred_data, K, max_week, season),
               ~ as.data.frame(fourier(..1$ILI, K = ..2,
                                       h = ..3 - 14 - 
                                         nrow(..1[..1$season == ..4, ])))),
        model == "Backfill" ~
          pmap(list(pred_data, K, max_week, season, backfill),
               ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                             h = ..3 - 14 - 
                                               nrow(..1[..1$season == ..4, ]))),
                       data.frame(backfill = ..5)) %>%
                 rename(`..1$backfill` = backfill)),
        model == "National Gtrends" ~
          pmap(list(pred_data, K, max_week, season, gtrend_forecast),
               ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                             h = ..3 - 14 - 
                                               nrow(..1[..1$season == ..4, ]))),
                       data.frame(hits = ..5)) %>%
                 rename(`..1$hits` = hits)),
        model == "Regional Gtrends" ~
          pmap(list(pred_data, K, max_week, season, reg_gtrend_forecast),
               ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                             h = ..3 - 14 - 
                                               nrow(..1[..1$season == ..4, ]))),
                       data.frame(region_hits = ..5)) %>%
                 rename(`..1$region_hits` = region_hits)),
        model == "FluVirus" ~
          pmap(list(pred_data, K, max_week, season, 
                    h1_per_forecast, h3_per_forecast, b_per_forecast),
               ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                             h = ..3 - 14 - 
                                               nrow(..1[..1$season == ..4, ]))),
                       data.frame(h1_per_samples = ..5,
                                  h3_per_samples = ..6, 
                                  b_per_samples = ..7)) %>%
                 rename(`..1$h1_per_samples` = h1_per_samples,
                        `..1$h3_per_samples` = h3_per_samples,
                        `..1$b_per_samples` = b_per_samples)),
        model == "National Gtrends & Backfill" ~
          pmap(list(pred_data, K, max_week, season, gtrend_forecast, backfill),
               ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                             h = ..3 - 14 - 
                                               nrow(..1[..1$season == ..4, ]))),
                       data.frame(hits = ..5,
                                  backfill = ..6)) %>%
                 rename(`..1$hits` = hits,
                        `..1$backfill` = backfill)),
        
        model == "Regional Gtrends & Backfill" ~
          pmap(list(pred_data, K, max_week, season, reg_gtrend_forecast, backfill),
               ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                             h = ..3 - 14 - 
                                               nrow(..1[..1$season == ..4, ]))),
                       data.frame(region_hits = ..5,
                                  backfill = ..6)) %>%
                 rename(`..1$region_hits` = region_hits,
                        `..1$backfill` = backfill)),
        model == "FluVirus & Backfill" ~
          pmap(list(pred_data, K, max_week, season, 
                    h1_per_forecast, h3_per_forecast, b_per_forecast,
                    backfill),
               ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                             h = ..3 - 14 - 
                                               nrow(..1[..1$season == ..4, ]))),
                       data.frame(h1_per_samples = ..5,
                                  h3_per_samples = ..6, 
                                  b_per_samples = ..7,
                                  backfill = ..8)) %>%
                 rename(`..1$h1_per_samples` = h1_per_samples,
                        `..1$h3_per_samples` = h3_per_samples,
                        `..1$b_per_samples` = b_per_samples,
                        `..1$backfill` = backfill)),
        model == "National Gtrends & FluVirus" ~
          pmap(list(pred_data, K, max_week, season, gtrend_forecast,
                    h1_per_forecast, h3_per_forecast, b_per_forecast),
               ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                             h = ..3 - 14 - 
                                               nrow(..1[..1$season == ..4, ]))),
                       data.frame(hits = ..5, 
                                  h1_per_samples = ..6,
                                  h3_per_samples = ..7, 
                                  b_per_samples = ..8)) %>%
                 rename(`..1$hits` = hits,
                        `..1$h1_per_samples` = h1_per_samples,
                        `..1$h3_per_samples` = h3_per_samples,
                        `..1$b_per_samples` = b_per_samples)),
        model == "Regional Gtrends & FluVirus" ~
          pmap(list(pred_data, K, max_week, season, reg_gtrend_forecast,
                    h1_per_forecast, h3_per_forecast, b_per_forecast),
               ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                             h = ..3 - 14 - 
                                               nrow(..1[..1$season == ..4, ]))),
                       data.frame(region_hits = ..5, 
                                  h1_per_samples = ..6,
                                  h3_per_samples = ..7, 
                                  b_per_samples = ..8)) %>%
                 rename(`..1$region_hits` = region_hits,
                        `..1$h1_per_samples` = h1_per_samples,
                        `..1$h3_per_samples` = h3_per_samples,
                        `..1$b_per_samples` = b_per_samples)),
        model == "National Gtrends, FluVirus, & Backfill" ~
          pmap(list(pred_data, K, max_week, season, gtrend_forecast,
                    h1_per_forecast, h3_per_forecast, b_per_forecast,
                    backfill),
               ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                             h = ..3 - 14 - 
                                               nrow(..1[..1$season == ..4, ]))),
                       data.frame(hits = ..5, 
                                  h1_per_samples = ..6,
                                  h3_per_samples = ..7, 
                                  b_per_samples = ..8,
                                  backfill = ..9)) %>%
                 rename(`..1$hits` = hits,
                        `..1$h1_per_samples` = h1_per_samples,
                        `..1$h3_per_samples` = h3_per_samples,
                        `..1$b_per_samples` = b_per_samples,
                        `..1$backfill` = backfill)),
        model == "Regional Gtrends, FluVirus, & Backfill" ~
          pmap(list(pred_data, K, max_week, season, reg_gtrend_forecast,
                    h1_per_forecast, h3_per_forecast, b_per_forecast,
                    backfill),
               ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                             h = ..3 - 14 - 
                                               nrow(..1[..1$season == ..4, ]))),
                       data.frame(region_hits = ..5, 
                                  h1_per_samples = ..6,
                                  h3_per_samples = ..7, 
                                  b_per_samples = ..8,
                                  backfill = ..9)) %>%
                 rename(`..1$region_hits` = region_hits,
                        `..1$h1_per_samples` = h1_per_samples,
                        `..1$h3_per_samples` = h3_per_samples,
                        `..1$b_per_samples` = b_per_samples,
                        `..1$backfill` = backfill))
      ),
      pred_results = pmap(
        list(pred_fit, pred_data, forecast_xreg, location, season, max_week),
        ~ fit_to_forecast(object = ..1,
                          xreg = ..3,
                          pred_data = ..2,
                          location = ..4,
                          season = ..5,
                          max_week = ..6,
                          npaths = 250)
      )
    ) %>%
    collect() %>%
    as.tibble()

  covar_scores <- covar_forecasts %>% ungroup() %>%
    select(season, model, location, week, pred_results) %>%
    mutate(pred_results = map2(pred_results, location,
                               ~ mutate(.x, location = .y) %>%
                                 normalize_probs())) %>%
    select(season, model, week, pred_results) %>%
    unnest() %>%
    nest(-season, -model, -week) %>%
    left_join(nested_truth, by = "season") %>%
    mutate(scores = map2(data, exp_truth,
                         ~ score_entry(.x, .y)),
           eval_scores = pmap(list(scores, eval_period, season),
                              ~ create_eval_scores(..1, ..2, ..3))) %>%
    select(season, model, eval_scores) %>%
    unnest() 
  
  save(covar_scores, file = paste0("Data/CV_covar_Scores_", 
                                  substr(this_season, 1, 4), ".Rdata"))
  
  rm(covar_forecasts, covar_model_parallel, cluster)
  
}
Sys.time() - start_time

load("Data/CV_covar_Scores_2010.Rdata")
covar_scores_1011 <- covar_scores
load("Data/CV_covar_Scores_2011.Rdata")
covar_scores_1112 <- covar_scores
load("Data/CV_covar_Scores_2012.Rdata")
covar_scores_1213 <- covar_scores
load("Data/CV_covar_Scores_2013.Rdata")
covar_scores_1314 <- covar_scores
load("Data/CV_covar_Scores_2014.Rdata")
covar_scores_1415 <- covar_scores
load("Data/CV_covar_Scores_2015.Rdata")
covar_scores_1516 <- covar_scores
load("Data/CV_covar_Scores_2016.Rdata")
covar_scores_1617 <- covar_scores
load("Data/CV_covar_Scores_2017.Rdata")
covar_scores_1718 <- covar_scores

# Determine best covariates for each region
best_covar_cv_scores <- bind_rows(covar_scores_1011, covar_scores_1112, covar_scores_1213,
                           covar_scores_1314, covar_scores_1415, covar_scores_1516,
                           covar_scores_1617, covar_scores_1718) %>%
  group_by(location, model) %>%
  summarize(avg_score = mean(score)) #%>%
  # group_by(location) %>%
  # filter(avg_score == max(avg_score)) %>%
  # ungroup() %>%
  # select(location, model)

# Manually select model for each region from results
# If multiple models similar (~ 0.01 diff), take model with most information
best_covar_cv <- tibble(location = c("HHS Region 1", "HHS Region 2", "HHS Region 3",
                                     "HHS Region 4", "HHS Region 5", "HHS Region 6",
                                     "HHS Region 7", "HHS Region 8", "HHS Region 9",
                                     "HHS Region 10", "US National"),
                        model = c("Regional Gtrends", "FluVirus & Backfill",
                                  "Regional Gtrends, FluVirus, & Backfill",
                                  "National Gtrends, FluVirus, & Backfill",
                                  "Backfill", "Regional Gtrends, FluVirus, & Backfill",
                                  "Regional Gtrends, FluVirus, & Backfill",
                                  "Regional Gtrends, FluVirus, & Backfill",
                                  "ARIMA only",
                                  "Regional Gtrends, FluVirus, & Backfill",
                                  "FluVirus & Backfill"))

save(best_covar_cv, file = "Data/CV_covar_terms.Rdata")


##### Create final fits and build forecasts #####
load("Data/CV_Transform_terms.Rdata")
load("data/CV_Fourier_terms.Rdata")
load("Data/CV_ARIMA_terms.Rdata")
load("Data/CV_covar_terms.Rdata")

# Set up model fits
load("Data/Final_CV_fits.Rdata")
# final_fits <- tibble(season = c("2010/2011", "2011/2012", "2012/2013",
#                                   "2013/2014", "2014/2015", "2015/2016",
#                                   "2016/2017", "2017/2018")) %>%
#   mutate(train_data = map(season,
#                           ~ filter(flu_data_merge, year <= as.numeric(substr(., 6, 9)),
#                                    season != paste0(substr(., 6, 9), "/",
#                                                     as.numeric(substr(., 6, 9)) + 1),
#                                    season != .))) %>%
#   unnest() %>%
#   # Nest by season and location
#   nest(-season, -location) %>%
#   # Create time series of ILI
#   mutate(data = map(data,
#                     ~ mutate(.x,
#                              ILI = ts(ILI, frequency = 52, start = c(2006, 40))))) %>%
#   # Merge best lambda, Fourier K value, and ARIMA structure by location
#   left_join(rename(best_transform_cv, transform = model), by = "location") %>%
#   mutate(lambda = unlist(map(data,
#                              ~ BoxCox.lambda(.$ILI))),
#          lambda = case_when(transform == "no_trans" ~ NA_real_,
#                             transform == "log" ~ 0,
#                             transform == "Box_Cox" ~ lambda)) %>%
#   left_join(best_k_cv, by = "location") %>%
#   left_join(best_arima_cv, by = "location") %>%
#   left_join(best_covar_cv, by = "location") %>%
#   # Fit models
#   mutate(fit = case_when(
#     model == "ARIMA only" ~
#       pmap(
#         list(data, arima_1, arima_2, arima_3, K, lambda),
#         ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
#                 xreg = fourier(..1$ILI, K = ..5),
#                 lambda = ..6)
#       ),
#     model == "National Gtrends" ~
#       pmap(
#         list(data, arima_1, arima_2, arima_3, K, lambda),
#         ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
#                 xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
#                              ..1$hits),
#                 lambda = ..6)
#       ),
#     model == "Regional Gtrends" ~
#       pmap(
#         list(data, arima_1, arima_2, arima_3, K, lambda),
#         ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
#                 xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
#                              ..1$region_hits),
#                 lambda = ..6)
#       ),
#     model == "FluVirus" ~
#       pmap(
#         list(data, arima_1, arima_2, arima_3, K, lambda),
#         ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
#                 xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
#                              ..1$h1_per_samples, ..1$h3_per_samples,
#                              ..1$b_per_samples),
#                 lambda = ..6)
#       ),
#     model == "Backfill" ~
#       pmap(
#         list(data, arima_1, arima_2, arima_3, K, lambda),
#         ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
#                 xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
#                              ..1$backfill),
#                 lambda = ..6)
#       ),
#     model == "FluVirus & Backfill" ~
#       pmap(
#         list(data, arima_1, arima_2, arima_3, K, lambda),
#         ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
#                 xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
#                              ..1$h1_per_samples, ..1$h3_per_samples,
#                              ..1$b_per_samples, ..1$backfill),
#                 lambda = ..6)
#       ),
#     model == "National Gtrends & FluVirus" ~
#       pmap(
#         list(data, arima_1, arima_2, arima_3, K, lambda),
#         ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
#                 xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
#                              ..1$hits, ..1$h1_per_samples,
#                              ..1$h3_per_samples, ..1$b_per_samples),
#                 lambda = ..6)
#       ),
#     model == "Regional Gtrends & FluVirus" ~
#       pmap(
#         list(data, arima_1, arima_2, arima_3, K, lambda),
#         ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
#                 xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
#                              ..1$region_hits, ..1$h1_per_samples,
#                              ..1$h3_per_samples, ..1$b_per_samples),
#                 lambda = ..6)
#       ),
#     model == "National Gtrends & Backfill" ~
#       pmap(
#         list(data, arima_1, arima_2, arima_3, K, lambda),
#         ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
#                 xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
#                              ..1$hits, ..1$backfill),
#                 lambda = ..6)
#       ),
#     model == "Regional Gtrends & Backfill" ~
#       pmap(
#         list(data, arima_1, arima_2, arima_3, K, lambda),
#         ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
#                 xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
#                              ..1$region_hits, ..1$backfill),
#                 lambda = ..6)
#       ),
#     model == "National Gtrends, FluVirus, & Backfill" ~
#       pmap(
#         list(data, arima_1, arima_2, arima_3, K, lambda),
#         ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
#                 xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
#                              ..1$hits, ..1$h1_per_samples,
#                              ..1$h3_per_samples, ..1$b_per_samples,
#                              ..1$backfill),
#                 lambda = ..6)
#       ),
#     model == "Regional Gtrends, FluVirus, & Backfill" ~
#       pmap(
#         list(data, arima_1, arima_2, arima_3, K, lambda),
#         ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
#                 xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
#                              ..1$region_hits, ..1$h1_per_samples,
#                              ..1$h3_per_samples, ..1$b_per_samples,
#                              ..1$backfill),
#                 lambda = ..6)
#       )
#   ))
# 
# save(final_fits, file = "Data/Final_CV_fits.Rdata")
  
final_forecast_data_setup <- crossing(season = c("2010/2011", "2011/2012", "2012/2013",
                                             "2013/2014", "2014/2015", "2015/2016",
                                             "2016/2017", "2017/2018"),
                                      week = c(40:73),
                                      location = unique(flu_data_merge$location)) %>%
  filter(week < 73 | season == "2014/2015") %>%
  mutate(
    epiweek = case_when(
      season == "2014/2015" & week > 53 ~ 
        as.numeric(paste0(substr(season, 6, 9), 
                          str_pad(week - 53, 2, "left", "0"))),
      week > 52 ~ 
        as.numeric(paste0(substr(season, 6, 9), 
                          str_pad(week - 52, 2, "left", "0"))),
      TRUE ~ as.numeric(paste0(substr(season, 1, 4), str_pad(week, 2, "left", "0")))
    ),
    prev_season = case_when(
      season == "2010/2011" ~ "2007/2008",
      TRUE ~ paste0(as.numeric(substr(season, 1, 4)) - 1, "/",
                    substr(season, 1, 4))
    )
  ) %>%
  inner_join(final_fits, 
             by = c("location", "season")) 

# Create and save forecast files
for(this_season in unique(final_forecast_data_setup$season)[-1]) {
  
  # Create folder if needed
  season_path <- paste0("Forecasts/", substr(this_season, 1, 4), "-",
                        substr(this_season, 6, 9), "/Dynamic Harmonic Model/")

  dir.create(season_path, showWarnings = FALSE)
  
  for(this_week in unique(final_forecast_data_setup$week[final_forecast_data_setup$season == this_season])) {

    temp <- final_forecast_data_setup %>%
      filter(season == this_season, week == this_week) %>%
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
        # Set up Fourier data for forecasting
        xreg = case_when(
          model == "ARIMA only" ~ 
            map2(pred_data, K,
                 ~ as.data.frame(fourier(.x$ILI, K = .y))),
          model == "National Gtrends" ~
            pmap(list(pred_data, K), 
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2)),
                         ..1$hits)),
          model == "Regional Gtrends" ~
            pmap(list(pred_data, K), 
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2)),
                         ..1$region_hits)),
          model == "FluVirus" ~ 
            pmap(list(pred_data, K),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2)),
                         ..1$h1_per_samples, ..1$h3_per_samples,
                         ..1$b_per_samples)),
          model == "Backfill" ~ 
            pmap(list(pred_data, K),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2)),
                         ..1$backfill)),
          model == "National Gtrends & Backfill" ~
            pmap(list(pred_data, K), 
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2)),
                         ..1$hits, ..1$backfill)),
          model == "Regional Gtrends & Backfill" ~
            pmap(list(pred_data, K), 
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2)),
                         ..1$region_hits, ..1$backfill)),
          model == "FluVirus & Backfill" ~ 
            pmap(list(pred_data, K),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2)),
                         ..1$h1_per_samples, ..1$h3_per_samples,
                         ..1$b_per_samples, ..1$backfill)),
          model == "National Gtrends & FluVirus" ~ 
            pmap(list(pred_data, K),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2)),
                         ..1$hits, ..1$h1_per_samples, 
                         ..1$h3_per_samples, ..1$b_per_samples)),
          model == "Regional Gtrends & FluVirus" ~ 
            pmap(list(pred_data, K),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2)),
                         ..1$region_hits, ..1$h1_per_samples, 
                         ..1$h3_per_samples, ..1$b_per_samples)),
          model == "National Gtrends, FluVirus, & Backfill" ~ 
            pmap(list(pred_data, K),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2)),
                         ..1$hits, ..1$h1_per_samples, 
                         ..1$h3_per_samples, ..1$b_per_samples,
                         ..1$backfill)),
          model == "Regional Gtrends, FluVirus, & Backfill" ~ 
            pmap(list(pred_data, K),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2)),
                         ..1$region_hits, ..1$h1_per_samples, 
                         ..1$h3_per_samples, ..1$b_per_samples,
                         ..1$backfill))
        ),
        max_week = ifelse(season == "2014/2015", 53, 52),
        # Fit models
        pred_fit = pmap(list(pred_data, fit, xreg),
                        ~ Arima(..1$ILI, xreg = ..3, model = ..2)),
        # Create data frame of xreg terms for forecasting
        gtrend_forecast = pmap(
          list(pred_data, location, season, prev_season, week, max_week),
          ~ tibble(hits = c(flu_data_merge %>%
                              filter(location == ..2, season == ..3, 
                                     order_week == ..5 + 1) %>%
                              pull(hits),
                            flu_data_merge %>%
                              filter(location == ..2, season == ..4, 
                                     order_week > ..5 + 1, 
                                     order_week < ..6 + 26) %>%
                              pull(hits) *
                              mean(flu_data_merge %>%
                                     filter(location == ..2, season == ..3,
                                            order_week <= ..5 + 1) %>%
                                     pull(hits) /
                                     flu_data_merge %>%
                                     filter(location == ..2, season == ..4,
                                            order_week <= ..5 + 1) %>%
                                     pull(hits))))
        ),
        reg_gtrend_forecast = pmap(
          list(pred_data, location, season, prev_season, week, max_week),
          ~ tibble(region_hits = c(flu_data_merge %>%
                                     filter(location == ..2, season == ..3, 
                                            order_week == ..5 + 1) %>%
                                     pull(region_hits),
                                   flu_data_merge %>%
                                     filter(location == ..2, season == ..4, 
                                            order_week > ..5 + 1, 
                                            order_week < ..6 + 26) %>%
                                     pull(region_hits) *
                                     mean(flu_data_merge %>%
                                            filter(location == ..2, season == ..3,
                                                   order_week <= ..5 + 1) %>%
                                            pull(region_hits) /
                                            flu_data_merge %>%
                                            filter(location == ..2, season == ..4,
                                                   order_week <= ..5 + 1) %>%
                                            pull(region_hits))))
        ),
        h1_per_forecast = pmap(
          list(pred_data, season, max_week),
          ~ rep(last(..1$h1_per_samples), ..3 - 14 - 
                  nrow(..1[..1$season == ..2, ]))
        ),
        h3_per_forecast = pmap(
          list(pred_data, season, max_week),
          ~ rep(last(..1$h3_per_samples), ..3 - 14 - 
                  nrow(..1[..1$season == ..2, ]))
        ),
        b_per_forecast = pmap(
          list(pred_data, season, max_week),
          ~ rep(last(..1$b_per_samples), ..3 - 14 - 
                  nrow(..1[..1$season == ..2, ]))
        ),
        backfill = pmap(
          list(data, week, max_week, season),
          function(data, week, max_week, season) {
            out <- numeric()
            for(i in week:(max_week + 24)) {
              
              if(season %in% c("2010/2011", "2011/2012", "2012/2013")) {
                temp <- data[data$order_week == i, ]
              } else{
                temp <- data[data$order_week == i &
                               !data$season %in% c("2010/2011", "2011/2012", "2012/2013"), ]
              }
              
              out <- c(out, predict(lm(backfill ~ year, 
                                       data = data[data$order_week == i, ]),
                                    data.frame(year = as.numeric(substr(season, 1, 4)))))
            }
            out
          } 
        ),
        forecast_xreg = case_when(
          model == "ARIMA only" ~ 
            pmap(list(pred_data, K, max_week, season),
                 ~ as.data.frame(fourier(..1$ILI, K = ..2,
                                         h = ..3 - 14 - 
                                           nrow(..1[..1$season == ..4, ])))),
          model == "Backfill" ~
            pmap(list(pred_data, K, max_week, season, backfill),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                               h = ..3 - 14 - 
                                                 nrow(..1[..1$season == ..4, ]))),
                         data.frame(backfill = ..5)) %>%
                   rename(`..1$backfill` = backfill)),
          model == "National Gtrends" ~
            pmap(list(pred_data, K, max_week, season, gtrend_forecast),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                               h = ..3 - 14 - 
                                                 nrow(..1[..1$season == ..4, ]))),
                         data.frame(hits = ..5)) %>%
                   rename(`..1$hits` = hits)),
          model == "Regional Gtrends" ~
            pmap(list(pred_data, K, max_week, season, reg_gtrend_forecast),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                               h = ..3 - 14 - 
                                                 nrow(..1[..1$season == ..4, ]))),
                         data.frame(region_hits = ..5)) %>%
                   rename(`..1$region_hits` = region_hits)),
          model == "FluVirus" ~
            pmap(list(pred_data, K, max_week, season, 
                      h1_per_forecast, h3_per_forecast, b_per_forecast),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                               h = ..3 - 14 - 
                                                 nrow(..1[..1$season == ..4, ]))),
                         data.frame(h1_per_samples = ..5,
                                    h3_per_samples = ..6, 
                                    b_per_samples = ..7)) %>%
                   rename(`..1$h1_per_samples` = h1_per_samples,
                          `..1$h3_per_samples` = h3_per_samples,
                          `..1$b_per_samples` = b_per_samples)),
          model == "National Gtrends & Backfill" ~
            pmap(list(pred_data, K, max_week, season, gtrend_forecast, backfill),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                               h = ..3 - 14 - 
                                                 nrow(..1[..1$season == ..4, ]))),
                         data.frame(hits = ..5,
                                    backfill = ..6)) %>%
                   rename(`..1$hits` = hits,
                          `..1$backfill` = backfill)),
          
          model == "Regional Gtrends & Backfill" ~
            pmap(list(pred_data, K, max_week, season, reg_gtrend_forecast, backfill),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                               h = ..3 - 14 - 
                                                 nrow(..1[..1$season == ..4, ]))),
                         data.frame(region_hits = ..5,
                                    backfill = ..6)) %>%
                   rename(`..1$region_hits` = region_hits,
                          `..1$backfill` = backfill)),
          model == "FluVirus & Backfill" ~
            pmap(list(pred_data, K, max_week, season, 
                      h1_per_forecast, h3_per_forecast, b_per_forecast,
                      backfill),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                               h = ..3 - 14 - 
                                                 nrow(..1[..1$season == ..4, ]))),
                         data.frame(h1_per_samples = ..5,
                                    h3_per_samples = ..6, 
                                    b_per_samples = ..7,
                                    backfill = ..8)) %>%
                   rename(`..1$h1_per_samples` = h1_per_samples,
                          `..1$h3_per_samples` = h3_per_samples,
                          `..1$b_per_samples` = b_per_samples,
                          `..1$backfill` = backfill)),
          model == "National Gtrends & FluVirus" ~
            pmap(list(pred_data, K, max_week, season, gtrend_forecast,
                      h1_per_forecast, h3_per_forecast, b_per_forecast),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                               h = ..3 - 14 - 
                                                 nrow(..1[..1$season == ..4, ]))),
                         data.frame(hits = ..5, 
                                    h1_per_samples = ..6,
                                    h3_per_samples = ..7, 
                                    b_per_samples = ..8)) %>%
                   rename(`..1$hits` = hits,
                          `..1$h1_per_samples` = h1_per_samples,
                          `..1$h3_per_samples` = h3_per_samples,
                          `..1$b_per_samples` = b_per_samples)),
          model == "Regional Gtrends & FluVirus" ~
            pmap(list(pred_data, K, max_week, season, reg_gtrend_forecast,
                      h1_per_forecast, h3_per_forecast, b_per_forecast),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                               h = ..3 - 14 - 
                                                 nrow(..1[..1$season == ..4, ]))),
                         data.frame(region_hits = ..5, 
                                    h1_per_samples = ..6,
                                    h3_per_samples = ..7, 
                                    b_per_samples = ..8)) %>%
                   rename(`..1$region_hits` = region_hits,
                          `..1$h1_per_samples` = h1_per_samples,
                          `..1$h3_per_samples` = h3_per_samples,
                          `..1$b_per_samples` = b_per_samples)),
          model == "National Gtrends, FluVirus, & Backfill" ~
            pmap(list(pred_data, K, max_week, season, gtrend_forecast,
                      h1_per_forecast, h3_per_forecast, b_per_forecast,
                      backfill),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                               h = ..3 - 14 - 
                                                 nrow(..1[..1$season == ..4, ]))),
                         data.frame(hits = ..5, 
                                    h1_per_samples = ..6,
                                    h3_per_samples = ..7, 
                                    b_per_samples = ..8,
                                    backfill = ..9)) %>%
                   rename(`..1$hits` = hits,
                          `..1$h1_per_samples` = h1_per_samples,
                          `..1$h3_per_samples` = h3_per_samples,
                          `..1$b_per_samples` = b_per_samples,
                          `..1$backfill` = backfill)),
          model == "Regional Gtrends, FluVirus, & Backfill" ~
            pmap(list(pred_data, K, max_week, season, reg_gtrend_forecast,
                      h1_per_forecast, h3_per_forecast, b_per_forecast,
                      backfill),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                               h = ..3 - 14 - 
                                                 nrow(..1[..1$season == ..4, ]))),
                         data.frame(region_hits = ..5, 
                                    h1_per_samples = ..6,
                                    h3_per_samples = ..7, 
                                    b_per_samples = ..8,
                                    backfill = ..9)) %>%
                   rename(`..1$region_hits` = region_hits,
                          `..1$h1_per_samples` = h1_per_samples,
                          `..1$h3_per_samples` = h3_per_samples,
                          `..1$b_per_samples` = b_per_samples,
                          `..1$backfill` = backfill))
        ),
        pred_results = pmap(
          list(pred_fit, pred_data, forecast_xreg, location, season, max_week),
          ~ fit_to_forecast(object = ..1,
                            xreg = ..3,
                            pred_data = ..2,
                            location = ..4,
                            season = ..5,
                            max_week = ..6,
                            npaths = 500)
        )
      ) %>%
      select(location, week, max_week, pred_results)
    
    EW <- case_when(temp$week[1] > temp$max_week[1] ~ 
                      str_pad(temp$week[1] - temp$max_week[1], 2, side = "left", 0),
                    TRUE ~ str_pad(temp$week[1], 2, "left", 0))
    
    temp %>%
      unnest() %>%
      select(-week, -max_week, -forecast_week) %>%
      bind_rows(generate_point_forecasts(.)) %>%
      select(location, target, type, unit, bin_start_incl, bin_end_notincl, value) %>%
      write_csv(path = paste0(season_path, "/EW", EW, ".csv"))
    
  }
  
}

