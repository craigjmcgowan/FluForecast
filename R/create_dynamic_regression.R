# Dynamic harmonic regression model 

library(tidyverse)
library(forecast)
library(lubridate)
library(FluSight)
library(MMWRweek)
library(parallel)
library(multidplyr)

# Load functions
source("R/utils.R")

# Load data
load("Data/ili.Rdata")
load("Data/virologic.Rdata")
load("Data/Gtrends.Rdata")

# Create clusters for use later in program
cluster <- create_cluster(cores = detectCores())

# Create truth for all seasons -------
load("Data/past_truth.Rdata")
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
# save(nested_truth, file = "Data/past_truth.Rdata")


# Combine datasets together -------
flu_data_merge <- select(ili_current, epiweek, ILI, year, week, season, location) %>%
  inner_join(select(virologic_combined, location, season, year, week, h1_per_samples,
                   h3_per_samples, b_per_samples),
            by = c("location", "season", "year", "week")) %>%
  inner_join(gtrend_US_flu_merge,
            by = c("season", "year", "week")) %>%
  # Remove 2008/2009 and 2009/2010 seasons due to pandemic activity
  filter(!season %in% c("2008/2009", "2009/2010")) %>%
  # Remove week 33 in 2014 so all seasons have 52 weeks - minimal activity
  filter(!(year == 2014 & week == 33)) %>%
  mutate(order_week = case_when(
    week < 40 & season == "2014/2015" ~ week + 53,
    week < 40 ~ week + 52,
    TRUE ~ week
  ))

### Decisions to make ###

## Number of Fourier terms
## Structure of ARIMA errors
## What virologic/Gtrend data to include
## If/how to include estimates of backfill\

##### Number of Fourier terms #####
# load("Data/fourier_fits.Rdata")
fourier_model_fits <- tibble(season = c("2010/2011", "2011/2012", "2012/2013", "2013/2014",
                                "2015/2016", "2016/2017", "2017/2018")) %>%
  mutate(train_data = map(season,
                          ~ filter(flu_data_merge, year <= as.numeric(substr(., 6, 9)),
                                   season != paste0(substr(., 6, 9), "/",
                                                    as.numeric(substr(., 6, 9)) + 1),
                                   season != .))) %>%
  unnest() %>%
  # Nest by season and location
  nest(-season, -location) %>%
  # Create time series of ILI
  mutate(data = map(data,
                   ~ mutate(.x,
                            ILI = ts(ILI, frequency = 52, start = c(2006, 40))))) %>%
  # Fit model
  mutate(fit_3 = map(data,
                     ~ auto.arima(.$ILI, xreg = fourier(.$ILI, K = 3),
                                  seasonal = FALSE, lambda = BoxCox.lambda(.$ILI))),
         fit_4 = map(data,
                     ~ auto.arima(.$ILI, xreg = fourier(.$ILI, K = 4),
                                  seasonal = FALSE, lambda = BoxCox.lambda(.$ILI))),
         fit_5 = map(data,
                     ~ auto.arima(.$ILI, xreg = fourier(.$ILI, K = 5),
                                  seasonal = FALSE, lambda = BoxCox.lambda(.$ILI))),
         fit_6 = map(data,
                     ~ auto.arima(.$ILI, xreg = fourier(.$ILI, K = 6),
                                  seasonal = FALSE, lambda = BoxCox.lambda(.$ILI))),
         fit_7 = map(data,
                     ~ auto.arima(.$ILI, xreg = fourier(.$ILI, K = 7),
                                  seasonal = FALSE, lambda = BoxCox.lambda(.$ILI))),
         fit_8 = map(data,
                     ~ auto.arima(.$ILI, xreg = fourier(.$ILI, K = 8),
                                  seasonal = FALSE, lambda = BoxCox.lambda(.$ILI))),
         fit_9 = map(data,
                     ~ auto.arima(.$ILI, xreg = fourier(.$ILI, K = 9),
                                  seasonal = FALSE, lambda = BoxCox.lambda(.$ILI))),
         fit_10 = map(data,
                     ~ auto.arima(.$ILI, xreg = fourier(.$ILI, K = 10),
                                  seasonal = FALSE, lambda = BoxCox.lambda(.$ILI))),
         fit_11 = map(data,
                     ~ auto.arima(.$ILI, xreg = fourier(.$ILI, K = 11),
                                  seasonal = FALSE, lambda = BoxCox.lambda(.$ILI)))) %>%
  select(-data) %>%
  gather(key = "model", value = "fit", fit_3:fit_11)

save(fourier_model_fits, file = "Data/fourier_fits.Rdata")

# Set up data for model fitting in parallel
fourier_model_data <- crossing(model = c("fit_3", "fit_4", "fit_5", "fit_6", "fit_7",
                                         "fit_8", "fit_9", "fit_10", "fit_11"),
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
  mutate(group = rep(1:length(cluster), length.out = nrow(.)))


# Set up party_df and load necessary libraries and functions 
fourier_by_group1 <- fourier_model_data %>%
  filter(season %in% c("2010/2011", "2011/2012", "2012/2013", "2013/2014")) %>%
  partition(group, cluster = cluster)

fourier_by_group1 %>% 
  cluster_library(c("tidyverse", "forecast", "lubridate", "FluSight", "MMWRweek")) %>% 
  cluster_assign_value("flu_data_merge", flu_data_merge) %>%
  cluster_assign_value("ili_init_pub_list", ili_init_pub_list) %>%
  cluster_assign_value("fit_to_forecast", fit_to_forecast) %>%
  cluster_assign_value("sample_predictive_trajectories_arima", 
                       sample_predictive_trajectories_arima)

# Create forecasts for Fourier terms
fourier_forecasts1 <- fourier_by_group1 %>%
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
Sys.time() - start_time
  
save(fourier_forecasts1, file = "Data/fourier_forecasts.Rdata") 
 
# Normalize probabilities and score forecasts 
fourier_scores <- fourier_forecasts %>% ungroup() %>%
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

save(fourier_scores, file = "Data/CV_Fourier_Scores.Rdata")

# Determine best K value for each region
best_k_cv <- fourier_scores %>%
  group_by(location, model) %>%
  summarize(avg_score = mean(score)) %>%
  filter(avg_score == min(avg_score)) %>%
  mutate(K = as.numeric(str_extract(model, "[1-9][0-9]|[0-9]"))) %>%
  ungroup() %>%
  select(location, K)

save(best_k_cv, file = "Data/CV_Fourier_terms.Rdata")

##### ARIMA structure for error terms #####

load("Data/CV_Fourier_terms.Rdata")
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
#   # Merge best Fourier K value by location
#   left_join(best_k_cv, by = "location") %>%
#   # Set up for parallel
#   mutate(group = rep(1:length(cluster), length.out = nrow(.)))
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
#     list(data, arima_1, arima_2, arima_3, K),
#     ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
#             xreg = fourier(..1$ILI, K = ..5),
#             lambda = BoxCox.lambda(..1$ILI))
#     )) %>%
#   select(-data) %>%
#   collect() %>%
#   as.tibble() %>%
#   ungroup()
# 
# save(arima_model_fits, file = "Data/arima_fits.Rdata")


# Set up data for forecast creation in parallel
arima_model_data_setup <- crossing(season = "2014/2015",
                                     # c("2010/2011", "2011/2012", "2012/2013",
                                     #          "2013/2014", "2014/2015", "2015/2016",
                                     #          "2016/2017", "2017/2018"),
                             arima_1 = 0:3,
                             arima_2 = 0:1,
                             arima_3 = 0:3,
                             week = 50, #c(43:71),
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
  left_join(select(arima_model_fits, -group), 
            by = c("location", "season", "arima_1", "arima_2", "arima_3")) %>%
  # Set up grouping for parallel
  mutate(group = rep(1:length(cluster), length.out = nrow(.)))
  

arima_model_data_parallel <- arima_model_data_setup %>%
  partition(group, cluster = cluster)

arima_model_data_parallel %>%
  cluster_library(c("tidyverse", "forecast", "lubridate", "FluSight", "MMWRweek")) %>% 
  cluster_assign_value("flu_data_merge", flu_data_merge) %>%
  cluster_assign_value("ili_init_pub_list", ili_init_pub_list) %>%
  cluster_assign_value("fit_to_forecast", fit_to_forecast) %>%
  cluster_assign_value("sample_predictive_trajectories_arima", 
                       sample_predictive_trajectories_arima)


# Create forecasts for different ARIMA structures
start_time <- Sys.time()
arima_forecasts <- arima_model_data_parallel %>%
  # filter(season == "2010/2011") %>%
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
                      npaths = 100))) %>%
  collect() %>%
  as.tibble() %>%
  ungroup()

Sys.time() - start_time

save(arima_forecasts, file = "Data/arima_forecasts.Rdata")

# Normalize probabilities and score forecasts 
arima_scores <- arima_forecasts %>%
  select(season, arima_1:arima_3, location, week, pred_results) %>%
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

save(arima_scores, file = "Data/CV_ARIMA_Scores.Rdata")

# Determine best ARIMA model for each region
best_arima_cv <- arima_scores %>%
  group_by(location, arima_1, arima_2, arima_3) %>%
  summarize(avg_score = mean(score)) %>%
  group_by(location) %>%
  filter(avg_score == min(avg_score)) %>%
  ungroup() %>%
  select(location, arima_1:arima_3)

save(best_arima_cv, file = "Data/CV_ARIMA_terms.Rdata")

##### Additional covariates #####
load("data/CV_Fourier_terms.Rdata")
load("Data/CV_ARIMA_terms.Rdata")

# Models to test:
# No covariates
# Google Trends data only
# Flu % only
# Gtrends and flu %

covar_model_fits <- crossing(season = "2017/2018",
                               # c("2010/2011", "2011/2012", "2012/2013",
                               #              "2013/2014", "2014/2015", "2015/2016",
                               #              "2016/2017", "2017/2018"),
                             model = c("ARIMA only", "Gtrends", "FluVirus", 
                                       "Gtrends & FluVirus")) %>%
  mutate(train_data = map(season,
                          ~ filter(flu_data_merge, year <= as.numeric(substr(., 6, 9)),
                                   season != paste0(substr(., 6, 9), "/",
                                                    as.numeric(substr(., 6, 9)) + 1),
                                   season != .))) %>%
  unnest() %>%
  # Nest by season and location
  nest(-season, -location, -model) %>%
  # Create time series of ILI
  mutate(data = map(data,
                   ~ mutate(.x,
                            ILI = ts(ILI, frequency = 52, start = c(2006, 40))))) %>%
  # Merge best Fourier K value and ARIMA structure by location
  left_join(best_k_cv, by = "location") %>%
  left_join(best_arima_cv, by = "location") %>%
  # Fit models
  mutate(fit = case_when(
    model == "ARIMA only" ~ 
      pmap(
        list(data, arima_1, arima_2, arima_3, K),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = fourier(..1$ILI, K = ..5),
                lambda = BoxCox.lambda(..1$ILI))
      ),
    model == "Gtrends" ~ 
      pmap(
        list(data, arima_1, arima_2, arima_3, K),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$hits),
                lambda = BoxCox.lambda(..1$ILI))
      ),
    model == "FluVirus" ~ 
      pmap(
        list(data, arima_1, arima_2, arima_3, K),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$h1_per_samples, ..1$h3_per_samples,
                             ..1$b_per_samples),
                lambda = BoxCox.lambda(..1$ILI))
      ),
    model == "Gtrends & FluVirus" ~
      pmap(
        list(data, arima_1, arima_2, arima_3, K),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$hits, ..1$h1_per_samples, 
                             ..1$h3_per_samples, ..1$b_per_samples),
                lambda = BoxCox.lambda(..1$ILI))
      )
    ))

# Set up data for forecast creation in parallel
covar_model_data_setup <- crossing(season = "2017/2018",
                                   # c("2010/2011", "2011/2012", "2012/2013",
                                   #          "2013/2014", "2014/2015", "2015/2016",
                                   #          "2016/2017", "2017/2018"),
                                   model = c("ARIMA only", "Gtrends", "FluVirus", 
                                             "Gtrends & FluVirus"),
                                   week = 50, #c(43:71),
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
  left_join(covar_model_fits, 
            by = c("location", "season", "model")) #%>%
  # Set up grouping for parallel
  # mutate(group = rep(1:length(cluster), length.out = nrow(.)))

# Set up dataset for forecasting in parallel
covar_model_parallel <- covar_model_data_setup %>%
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
      model == "Gtrends" ~
        pmap(list(pred_data, K), 
             ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2)),
                     ..1$hits)),
      model == "FluVirus" ~ 
        pmap(list(pred_data, K),
             ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2)),
                     ..1$h1_per_samples, ..1$h3_per_samples,
                     ..1$b_per_samples)),
      model == "Gtrends & FluVirus" ~ 
        pmap(list(pred_data, K),
             ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2)),
                   ..1$hits, ..1$h1_per_samples, 
                   ..1$h3_per_samples, ..1$b_per_samples))
      ),
    max_week = ifelse(season == "2014/2015", 53, 52),
    # Fit models
    pred_fit = pmap(list(pred_data, fit, xreg),
                    ~ Arima(..1$ILI, xreg = ..3, model = ..2)),
    # Create data frame of xreg terms for forecasting
    gtrend_forecast = pmap(
      list(pred_data, season, week, max_week),
      ~ tibble(hits = c(flu_data_merge %>%
                          filter(season == ..2, order_week == ..3 + 1) %>%
                          slice(1) %>%
                          pull(hits),
                        snaive(ts(..1$hits, 
                                  frequency = 52, 
                                  start = c(2004, 40)))$mean[2:(..4 - 17 - 
                                                                  nrow(..1[..1$season == ..2, ]))] *
                          flu_data_merge %>%
                            filter(season == ..2, order_week == ..3 + 1) %>%
                            slice(1) %>%
                            pull(hits) / 
                          snaive(ts(..1$hits, 
                                    frequency = 52, 
                                    start = c(2004, 40)))$mean[1]))
    ),
    h1_per_forecast = pmap(
      list(pred_data, season, max_week),
      ~ rep(last(..1$h1_per_samples), ..3 - 17 - 
              nrow(..1[..1$season == ..2, ]))
    ),
    h3_per_forecast = pmap(
      list(pred_data, season, max_week),
      ~ rep(last(..1$h3_per_samples), ..3 - 17 - 
              nrow(..1[..1$season == ..2, ]))
    ),
    b_per_forecast = pmap(
      list(pred_data, season, max_week),
      ~ rep(last(..1$b_per_samples), ..3 - 17 - 
              nrow(..1[..1$season == ..2, ]))
    ),
    forecast_xreg = case_when(
      model == "ARIMA only" ~ 
        pmap(list(pred_data, K, max_week, season),
             ~ as.data.frame(fourier(..1$ILI, K = ..2,
                                     h = ..3 - 17 - 
                                       nrow(..1[..1$season == ..4, ])))),
      model == "Gtrends" ~
        pmap(list(pred_data, K, max_week, season, gtrend_forecast),
             ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                           h = ..3 - 17 - 
                                             nrow(..1[..1$season == ..4, ]))),
                     data.frame(hits = ..5)) %>%
               rename(`..1$hits` = hits)),
      model == "FluVirus" ~
        pmap(list(pred_data, K, max_week, season, 
                  h1_per_forecast, h3_per_forecast, b_per_forecast),
             ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                           h = ..3 - 17 - 
                                             nrow(..1[..1$season == ..4, ]))),
                     data.frame(h1_per_samples = ..5,
                                h3_per_samples = ..6, 
                                b_per_samples = ..7)) %>%
               rename(`..1$h1_per_samples` = h1_per_samples,
                      `..1$h3_per_samples` = h3_per_samples,
                      `..1$b_per_samples` = b_per_samples)),
      model == "Gtrends & FluVirus" ~
        pmap(list(pred_data, K, max_week, season, gtrend_forecast,
                  h1_per_forecast, h3_per_forecast, b_per_forecast),
             ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                           h = ..3 - 17 - 
                                             nrow(..1[..1$season == ..4, ]))),
                     data.frame(hits = ..5, 
                                h1_per_samples = ..6,
                                h3_per_samples = ..7, 
                                b_per_samples = ..8)) %>%
               rename(`..1$hits` = hits,
                      `..1$h1_per_samples` = h1_per_samples,
                      `..1$h3_per_samples` = h3_per_samples,
                      `..1$b_per_samples` = b_per_samples))
    ),
    pred_results = pmap(
      list(pred_fit, pred_data, forecast_xreg, location, season, max_week),
      ~ fit_to_forecast(object = ..1,
                        xreg = ..3,
                        pred_data = ..2,
                        location = ..4,
                        season = ..5,
                        max_week = ..6,
                        npaths = 20)
    )
  )# %>%
  # collect() %>%
  # as.tibble()


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

save(arima_scores, file = "Data/CV_ARIMA_Scores.Rdata")

# Determine best ARIMA model for each region
best_covar_cv <- covar_scores %>%
  group_by(location, model) %>%
  summarize(avg_score = mean(score)) %>%
  group_by(location) %>%
  filter(avg_score == min(avg_score)) %>%
  ungroup() %>%
  select(location, model)

save(best_arima_cv, file = "Data/CV_ARIMA_terms.Rdata")

