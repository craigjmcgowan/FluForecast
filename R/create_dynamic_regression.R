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

#### Number of Fourier terms #####
load("Data/fourier_fits.Rdata")
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
#   # Create time series of ILI
#   mutate(data = map(data,
#                    ~ mutate(.x,
#                             ILI = ts(ILI, frequency = 52, start = c(2006, 40))))) %>%
#   # Fit model
#   mutate(fit_5 = map(data,
#                    ~ auto.arima(.$ILI, xreg = fourier(.$ILI, K = 5),
#                                 seasonal = FALSE, lambda = -BoxCox.lambda(.$ILI))),
#          fit_6 = map(data,
#                      ~ auto.arima(.$ILI, xreg = fourier(.$ILI, K = 6),
#                                   seasonal = FALSE, lambda = -BoxCox.lambda(.$ILI))),
#          fit_7 = map(data,
#                      ~ auto.arima(.$ILI, xreg = fourier(.$ILI, K = 7),
#                                   seasonal = FALSE, lambda = -BoxCox.lambda(.$ILI))),
#          fit_8 = map(data,
#                      ~ auto.arima(.$ILI, xreg = fourier(.$ILI, K = 8),
#                                   seasonal = FALSE, lambda = -BoxCox.lambda(.$ILI))),
#          fit_9 = map(data,
#                      ~ auto.arima(.$ILI, xreg = fourier(.$ILI, K = 9),
#                                   seasonal = FALSE, lambda = -BoxCox.lambda(.$ILI)))) %>%
#   select(-data) %>%
#   gather(key = "model", value = "fit", fit_5:fit_9)
# 
# save(fourier_model_fits, file = "Data/fourier_fits.Rdata")

# Set up data for model fitting in parallel
fourier_model_data <- crossing(model = c("fit_5", "fit_6", "fit_7", "fit_8", "fit_9"),
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
fourier_by_group <- fourier_model_data %>%
  partition(group, cluster = cluster)

fourier_by_group %>% 
  cluster_library(c("tidyverse", "forecast", "lubridate", "FluSight", "MMWRweek")) %>% 
  cluster_assign_value("flu_data_merge", flu_data_merge) %>%
  cluster_assign_value("ili_init_pub_list", ili_init_pub_list) %>%
  cluster_assign_value("fit_to_forecast", fit_to_forecast) %>%
  cluster_assign_value("sample_predictive_trajectories_arima", 
                       sample_predictive_trajectories_arima)



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
  mutate(K = as.numeric(str_extract(model, "[0-9]"))) %>%
  ungroup() %>%
  select(location, K)

save(best_k_cv, file = "Data/CV_Fourier_terms.Rdata")

### ARIMA structure for error terms ####

load("Data/CV_Fourier_terms.Rdata")
# load("Data/arima_fits.Rdata")

arima_model_fit_data <- crossing(season = "2010/2011",
                                #    c("2010/2011", "2011/2012", "2012/2013", "2013/2014",
                                # "2015/2016", "2016/2017", "2017/2018"),
                                arima_1 = 0:3,
                                arima_2 = 0:1,
                                arima_3 = 0:3) %>%
  mutate(train_data = map(season,
                          ~ filter(flu_data_merge, year <= as.numeric(substr(., 6, 9)),
                                   season != paste0(substr(., 6, 9), "/",
                                                    as.numeric(substr(., 6, 9)) + 1),
                                   season != .))) %>%
  unnest() %>%
  # Nest by season and location
  nest(-season, -location, -arima_1, -arima_2, -arima_3) %>%
  # Create time series of ILI
  mutate(data = map(data,
                   ~ mutate(.x,
                            ILI = ts(ILI, frequency = 52, start = c(2006, 40))))) %>%
  # Merge best Fourier K value by location
  left_join(best_k_cv, by = "location") %>%
  # Set up for parallel
  mutate(group = rep(1:length(cluster), length.out = nrow(.)))

# Set up clusters
arima_fit_parallel <- arima_model_fit_data %>%
  partition(group, cluster = cluster)
arima_fit_parallel %>%
  cluster_library(c("tidyverse", "forecast", "lubridate", "FluSight", "MMWRweek")) 

# Fit ARIMA models
arima_model_fits <- arima_fit_parallel %>%
  mutate(fit = pmap(
    list(data, arima_1, arima_2, arima_3, K),
    ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
            xreg = fourier(..1$ILI, K = ..5),
            lambda = -BoxCox.lambda(..1$ILI))
    )) %>%
  select(-data) %>%
  collect() %>%
  as.tibble() %>%
  ungroup()

save(arima_model_fits, file = "Data/arima_fits.Rdata")


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
arima_forecasts_1011 <- arima_model_data_parallel %>%
  filter(season == "2010/2011") %>%
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

### Additional covariates ###
load("data/CV_Fourier_terms.Rdata")
load("Data/CV_ARIMA_terms.Rdata")


# snaive_gtrend <- snaive(flu_data_merge$ILI)
# autoplot(snaive_gtrend)
# checkresiduals(snaive_gtrend)
# gglagplot(flu_data_merge$ILI)
# ggsubseriesplot(flu_data_merge$ILI)
# ggseasonplot(gtrend_US_ts)
# ?auto.arima
# ggAcf(flu_data_merge$ILI)
# diff(diff(flu_data_merge$ILI, 52), 1) %>% autoplot()
# mstl(ili_US_ts) %>%
#   forecast(method = "arima", bootstrap = TRUE) %>%
#   autoplot()
# ma2x52 <- ma(ili_US_ts, order = 52, centre = TRUE)
# autoplot(ili_US_ts) +
#   autolayer(ma2x52)
# decompose(flu_data_merge$ILI) %>%
#   autoplot()
# seasonal::seas(ili_US_ts, x11 = "")
# autoplot(gtrend_ILI_no2009$ILI)
# 
# BoxCox.lambda(ili_US_ts)
# 
# plots <- list()
# for (i in 1:6) {
#   fit <- auto.arima(flu_data_merge$ILI,
#                       xreg = fourier(flu_data_merge$ILI, K = i),
#                     seasonal = FALSE,
#                     lambda = -1)
#   plots[[paste(i)]] <- autoplot(forecast(fit,
#                                          xreg = )) +
#     xlab(paste("k=",i,"   AICC=",round(fit[["aicc"]],2))) +
#     ylab("")
# }
# 
# forecast(fit, xreg = fourier(gtrend_observed$ILI, K = 1, h = 30)) %>%
#   autoplot()
# 
# checkresiduals(fit)
# gridExtra::grid.arrange(
#   plots[[1]],plots[[2]],plots[[3]],
#   plots[[4]],plots[[5]],plots[[6]],
#   nrow=3)
# 
# ggplot(data = gtrend_ILI_no2009, aes(x = hits, y = ILI)) +
#   geom_point()
# ?simulate.Arima
# 
# 
# fit <- auto.arima(gtrend_observed$ILI,
#                   xreg = cbind(fourier(gtrend_observed$ILI, K = 4),
#                                gtrend_observed$hits),
#                   seasonal = FALSE,
#                   lambda = -1)
# 
# test <- sample_predictive_trajectories_arima(
#   fit,
#   h = 52,
#   xreg = cbind(fourier(gtrend_observed$ILI, K = 4, h = 52),
#                snaive(ts(gtrend_observed$hits, frequency = 52, start = c(2006, 40)), h = 52)$mean)
# )
# ?snaive()
# hist(test[, 4])
# dim(test)
# (fit <- auto.arima(gtrend_observed$ILI, lambda = 0
#                    xreg = cbind(fourier(ili_US_ts, 3),
#                                 gtrend_ILI_merge[1:618, c("hits", "h1_per_samples",
#                                                        "h3_per_samples")]),
#                    lambda = -1))
# (fit)
# 
# checkresiduals(fit)
# 
# BoxCox(ili_US_ts, -1) %>% autoplot()
# 
# forecast(fit, 
#          xreg = cbind(fourier(ili_US_ts, 3, 24),
#                            snaive(ts(gtrend_ILI_merge$hits[1:723], frequency = 52, start = c(2004, 40)))$mean[1:24],
#                            snaive(ts(gtrend_ILI_merge$h1_per_samples[1:723], frequency = 52, start = c(2004, 40)))$mean[1:24],
#                            snaive(ts(gtrend_ILI_merge$h3_per_samples[1:723], frequency = 52, start = c(2004, 40)))$mean[1:24]),
#          bootstrap = TRUE) %>%
#   autoplot()
# 
# forecast(fit, level = 80) %>% autoplot()
# 
# 
# ## For future Google Trends - compare the 1 wk known value to what was observed at the same
# # week last year. Adjust the 2, 3, 4 seasonal naive values by the ratio of the two 1 wk values