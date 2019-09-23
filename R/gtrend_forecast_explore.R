### Explore different methods of predicting Gtrend for next three weeks in models

library(tidyverse)
library(forecast)
library(zoo)

source('R/utils.R')

gtrend <- readRDS('Data/gtrend.Rds') %>%
  mutate(location = case_when(
    location == "Massachusetts" ~ "HHS Region 1",
    location == "New York" ~ "HHS Region 2",
    location == "Pennsylvania" ~ "HHS Region 3",
    location == "Florida" ~ "HHS Region 4",
    location == "Illinois" ~ "HHS Region 5",
    location == "Texas" ~ "HHS Region 6",
    location == "Missouri" ~ "HHS Region 7",
    location == "Colorado" ~ "HHS Region 8",
    location == "California" ~ "HHS Region 9",
    location == "Washington" ~ "HHS Region 10",
    TRUE ~ location
  )) %>%
  filter(!location %in% state.name, !season %in% c("2008/2009", "2009/2010"),
         !(year == 2015 & week == 33)) %>%
  # Create rolling 52 week mean for comparison
  group_by(location, season) %>%
  arrange(year, week, .by_group = TRUE) %>%
  mutate(year_avg = cummean(hits)) %>%
  ungroup() %>%
  # Create ratios and order_week variable
  mutate(ratio = year_avg / lag(year_avg, 52),
         order_week = week_inorder(week, season))




# model_fits <- tibble(season = c("2010/2011", "2011/2012", "2012/2013",
#                                 "2013/2014", "2014/2015", "2015/2016",
#                                 "2016/2017", "2017/2018")) %>%
#   mutate(train_data = map(season,
#                           ~ filter(gtrend, year <= as.numeric(substr(., 6, 9)),
#                                    season != paste0(substr(., 6, 9), "/",
#                                                     as.numeric(substr(., 6, 9)) + 1),
#                                    season != .) %>%
#                             select(-season, -year_avg, -ratio))) %>%
#   unnest(cols = c(train_data)) %>%
#   # Nest by season and location
#   nest(data = c(date, week, year, hits)) %>%
#   # Create time series of ILI and numerical value of lambda
#   mutate(data = map(data,
#                     ~ mutate(.x,
#                              hits = ts(hits, frequency = 52, start = c(2006, 40))))) %>%
#   # Fit models
#   mutate(
#     arima_model = map(data, ~ auto.arima(.$hits)),
#     snaive_model = map(data, ~ snaive(.$hits))
#   ) %>%
#   select(-data)
# 
# saveRDS(model_fits, "Data/gtrend_explore_fits.RDS")

model_fits <- readRDS("Data/gtrend_explore_fits.RDS")

model_forecasts <- crossing(season = c("2010/2011", "2011/2012", "2012/2013", 
                                       "2013/2014", "2014/2015", "2015/2016", 
                                       "2016/2017", "2017/2018"),
                            week = c(43:70),
                            location = unique(gtrend$location)) %>%
  # Set prediction data
  mutate(pred_data = pmap(select(., location, season, week),
                         ~ filter(gtrend, location == ..1,
                                  year <= as.numeric(substr(..2, 6, 9)),
                                  season != paste0(substr(..2, 6, 9), "/",
                                                   as.numeric(substr(..2, 6, 9)) + 1),
                                  season != ..2 | order_week %in% 40:..3) %>%
                           mutate(hits = ts(hits, frequency = 52, 
                                            start = c(2006, 40))) %>%
                           select(-season, -location, -year_avg)),
         truth = pmap(select(., location, season, week),
                      ~ filter(gtrend, location == ..1,
                               season == ..2,
                               order_week > ..3,
                               order_week <= 74))) %>%
  # Join model fit objects 
  left_join(model_fits, by = c("season", "location")) %>%
  # Create predictions
  mutate(snaive_forecast = pmap(select(., snaive_model, week, pred_data, truth),
                                ~ tibble(pred_week = (..2 + 1):74,
                                         snaive = ..1$mean[(..2 - 39):34] * tail(..3$ratio, 1),
                                         truth = ..4$hits,
                                         AE = truth - snaive)),
         arima_forecast = pmap(select(., arima_model, week, pred_data, truth),
                               ~ tibble(pred_week = (..2 + 1):74,
                                        arima = forecast(..3$hits, model = ..1,
                                                         h = (74 - ..2))$mean,
                                        truth = ..4$hits,
                                        AE = truth - arima))) %>%
  # Drop model fits
  select(-arima_model, -snaive_model)
  
# Evaluate forecasts - use MAE since point values are what matters for forecast
model_scores <- model_forecasts %>%
  mutate(snaive_MAE = map_dbl(snaive_forecast,
                              ~ mean(abs(.$AE))),
         arima_MAE = map_dbl(arima_forecast,
                             ~ mean(abs(.$AE)))) %>%
  select(location, season, week, snaive_MAE, arima_MAE)

# Average across different groupings
model_scores %>%
  summarize(snaive = mean(snaive_MAE),
            arima = mean(arima_MAE))

model_scores %>%
  group_by(location) %>%
  summarize(snaive = mean(snaive_MAE),
            arima = mean(arima_MAE))

model_scores %>%
  group_by(season) %>%
  summarize(snaive = mean(snaive_MAE),
            arima = mean(arima_MAE))

week_avg<-model_scores %>%
  group_by(week) %>%
  summarize(snaive = mean(snaive_MAE),
            arima = mean(arima_MAE))
