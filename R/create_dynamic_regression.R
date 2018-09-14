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
# load("Data/temp_forecasts.Rdata")

# Combine datasets together
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
## If/how to include estimates of backfill

## Evaluate based on CDC log score


## Plan of action:

## Set up model of some kind to produce CDC-style forecast sequentially for each week of the season
## Set up way to score it
## Scale up to evaluate lots of different model options


# Predictions for 2017-2018 season - simple model, just Fourier terms, no covariates
# Train model based on prior data
model_fits <- filter(flu_data_merge, year < 2018, season != "2017/2018") %>%
  # Nest by location
  nest(-location) %>%
  # Create time series of ILI
  mutate(data = map(data,
                   ~ mutate(.x,
                            ILI = ts(ILI, frequency = 52, start = c(2006, 40))))) %>%
  # Fit model
  mutate(fit = map(data,
                   ~ auto.arima(.$ILI, xreg = fourier(.$ILI, K = 8),
                                seasonal = FALSE, lambda = -BoxCox.lambda(.$ILI))),
         model = "Test Model",
         season = "2017/2018")

forecasts <- list()

this_season <- "2017/2018"
this_model <- "Test Model"

# Set up data for model fitting in parallel
model_data <- expand.grid(model = c("Test Model"),
                          season = c("2017/2018"),
                          week = c(43:70),
                          location = unique(flu_data_merge$location),
                          stringsAsFactors = FALSE) %>%
  mutate(epiweek = case_when(
      season == "2014/2015" & week > 53 ~ 
        as.numeric(paste0(substr(season, 6, 9), 
                          str_pad(week - 53, 2, "left", "0"))),
      week > 52 ~ 
        as.numeric(paste0(substr(season, 6, 9), 
                          str_pad(week - 52, 2, "left", "0"))),
      TRUE ~ as.numeric(paste0(substr(season, 1, 4), str_pad(week, 2, "left", "0")))
    ),
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
                       mutate(ILI = ts(ILI, frequency = 52, start = c(2006, 40))))) %>%
  left_join(select(model_fits, -data), by = c("season", "location", "model")) %>%
  mutate(xreg = map(pred_data, ~ fourier(.$ILI, K = 8))) %>%
  # Set up grouping for parallel
  mutate(group = rep(1:4, length.out = nrow(.)))

cluster <- create_cluster(cores = detectCores())

by_group <- model_data %>%
  partition(group, cluster = cluster)

by_group %>%
  cluster_library(c("tidyverse", "forecast", "lubridate", "FluSight", "MMWRweek")) %>%
  cluster_assign_value("fit_to_forecast", fit_to_forecast) %>%
  cluster_assign_value("sample_predictive_trajectories_arima", 
                       sample_predictive_trajectories_arima)
  
# Run code to create forecasts in parallel
start_time <- Sys.time()
forecasts <- by_group %>%
  # Fit ARIMA model based on most recent data using prior fit
  mutate(pred_fit = pmap(list(pred_data, fit, xreg),
                        ~ Arima(..1$ILI, xreg = ..3, model = ..2))) %>%
  # Create predicted results
  mutate(pred_results = pmap(
    list(pred_fit, pred_data, location, season),
    ~ fit_to_forecast(object = ..1,
                      xreg = fourier(..2$ILI, K = 8,
                                     h = 35 - nrow(filter(..2, season == ..4))),
                      pred_data = ..2,
                      location = ..3,
                      max_week = 52,
                      npaths = 500))) %>%
  collect() %>%
  as.tibble() %>%
  ungroup()

Sys.time() - start_time


forecasts2 <- forecasts %>% ungroup() %>%
  select(season, model, location, week, pred_results) %>%
  unnest()
  mutate(pred_results = map2(pred_results, location,
                             ~ mutate(.x, location = .y) %>%
                               normalize_probs())) %>%
  select(season, model, location, week, pred_results) 

test1 <- model_data$pred_data[[1]]
test2 <- model_data$pred_data[[2]]

  select(location, pred_results) %>%
  unnest() %>%
  normalize_probs()
names(forecasts2)
test <- forecasts2$test_results[[1]]


Sys.time()-start_time
save(forecasts, forecasts_list, file = "Data/temp_forecasts.Rdata")


# Create truth to score models against
nested_truth <- ili_current %>%
  filter(year >= 2010, season != "2009/2010") %>%
  select(season, location, week, ILI) %>%
  nest(-season) %>%
  mutate(truth = map2(season, data,
                      ~ create_truth(fluview = FALSE, year = substr(.x, 1, 4),
                                     weekILI = .y)),
         eval_period = pmap(list(data, truth, season),
                             ~ create_eval_period(..1, ..2, ..3)),
         exp_truth = map(truth,
                         ~ expand_truth(.))) %>%
  select(-data)

# Create nested df of forecasts and score
forecasts_df <- bind_rows(map(modify_depth(forecasts, 2, bind_rows), 
                              bind_rows, .id = "team"),
                          .id = "season") %>%
  mutate(forecast_week = as.numeric(forecast_week)) %>%
  mutate(extra_week = forecast_week) %>%
  nest(-season, -team, -extra_week) %>%
  left_join(nested_truth, by = "season") %>%
  mutate(scores = map2(data, exp_truth,
                       ~ score_entry(.x, .y)),
         eval_scores = pmap(list(scores, eval_period, season),
                            ~ create_eval_scores(..1, ..2, ..3))) 
  


scores <- forecasts_df %>%
  select(season, team, eval_scores) %>%
  unnest() %>%
  group_by(season, team) %>%
  summarize(mean_score = mean(score),
            skill = exp(mean_score))




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