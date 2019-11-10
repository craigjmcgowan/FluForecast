# Dynamic harmonic regression model 


# Create clusters for use later in program
library(tidyverse)
library(forecast)
library(lubridate)
library(FluSight)
library(multidplyr)
library(MMWRweek)

# Load functions
source("R/utils.R")

# Load data
# ili_backfill <- readRDS('Data/ili_backfill.RDS') %>%
#   filter(location %in% state.name)
# ili_backfill_month_avg <- readRDS('Data/ili_backfill_month_avg.RDS') %>%
#   filter(location %in% state.name)
ili_backfill_densities <- readRDS('Data/ili_backfill_densities.Rds') %>%
  filter(location %in% state.name) %>%
  select(-season)
ili_current <- readRDS('Data/ili_current.RDS') %>%
  filter(location %in% state.name)
ili_init_pub_list <- readRDS('Data/ili_init_pub_list.RDS') %>%
  lapply(function(x) dplyr::filter(x, location %in% state.name))
virologic_combined <- readRDS('Data/virologic.RDS') %>%
  filter(location %in% state.name)
gtrend<- readRDS('Data/gtrend.RDS') %>%
  filter(location %in% c(state.name, "US National"),
         !season %in% c("2008/2009", "2009/2010"),
         !(year == 2015 & week == 33))

# Create clusters for use later in program
cluster <- create_cluster(cores = parallel::detectCores() - 1)

# Create truth for all seasons and combine datasets -------
ILI_1415 <- ili_init_pub_list[["201740"]] %>%
  filter(week >= 40 | week < 25, season == "2014/2015",
         location %in% state.name) %>%
  mutate(week = week_inorder(week, season)) %>%
  arrange(location, week) %>%
  mutate(week = week_reset(week, season)) %>%
  select(season, location, week, ILI)

ILI_1516 <- ili_init_pub_list[["201740"]] %>%
  filter(week >= 40 | week < 25, season == "2015/2016",
         location %in% state.name) %>%
  mutate(week = week_inorder(week, season)) %>%
  arrange(location, week) %>%
  mutate(week = week_reset(week, season)) %>%
  select(season, location, week, ILI)

ILI_1617 <-  ili_init_pub_list[["201740"]] %>%
  filter(week >= 40 | week < 25, season == "2016/2017",
         location %in% state.name) %>%
  mutate(week = week_inorder(week, season)) %>%
  arrange(location, week) %>%
  mutate(week = week_reset(week, season)) %>%
  select(season, location, week, ILI)

ILI_1718 <-  ili_init_pub_list[["201828"]] %>%
  filter(week >= 40 | week < 25, season == "2017/2018",
         location %in% state.name) %>%
  mutate(week = week_inorder(week, season)) %>%
  arrange(location, week) %>%
  mutate(week = week_reset(week, season)) %>%
  select(season, location, week, ILI)

ILI_1819 <-  ili_init_pub_list[["201928"]] %>%
  filter(week >= 40 | week < 25, season == "2018/2019",
         location %in% state.name) %>%
  mutate(week = week_inorder(week, season)) %>%
  arrange(location, week) %>%
  mutate(week = week_reset(week, season)) %>%
  select(season, location, week, ILI)

state_nested_truth <- bind_rows(ILI_1415, ILI_1516, ILI_1617, 
                                ILI_1718, ILI_1819) %>%
  nest(data = c(location, week, ILI)) %>%
  mutate(include_53 = ifelse(season %in% c("2014/2015"), TRUE, FALSE),
         truth = map2(season, data,
                      ~ create_truth(fluview = FALSE, year = as.numeric(substr(.x, 1, 4)),
                                     weekILI = .y, challenge = "state_ili"))) %>%
  select(-data)


state_flu_data_merge <- select(ili_current, epiweek, ILI, year, week, season, location) %>%
  # inner_join(select(virologic_combined, location, season, year, week, cum_h1per,
  #                  cum_h3per, cum_bper),
  #           by = c("location", "season", "year", "week")) %>%
  # Add Google Trends data by location
  right_join(gtrend %>%
               filter(year > 2010 | (year == 2010 & week >= 40),
                      location %in% state.name) %>%
               filter(location != "Florida" &
                        (location != "Louisiana" | season %in% c('2016/2017', '2017/2018', '2018/2019'))) %>%
               rename(region_hits = hits) %>%
               select(-date),
            by = c("location", "season", "year", "week")) %>%
  right_join(filter(gtrend, year > 2010 | (year == 2010 & week >= 40),
                    location == "US National") %>%
               select(-location, -date),
            by = c("season", "year", "week")) %>%
  # Remove week 33 in 2015 so all seasons have 52 weeks - minimal activity
  filter(!(year == 2015 & week == 33)) %>%
  mutate(order_week = week_inorder(week, season),
         month = month(MMWRweek2Date(year, week))) 

### Decisions to make ###

## Number of Fourier terms
## Structure of ARIMA errors
## What virologic/Gtrend data to include
## If/how to include estimates of backfill\

##### Number of Fourier terms #####
fourier_model_fits <- readRDS("Data/state_fourier_fits.Rds")

fourier_model_fit_data <- tibble(season = c("2014/2015", "2015/2016", "2016/2017", 
                                            "2017/2018", "2018/2019")) %>%
  mutate(train_data = map(season,
                          ~ filter(state_flu_data_merge, year <= as.numeric(substr(., 6, 9)),
                                   season != paste0(substr(., 6, 9), "/",
                                                    as.numeric(substr(., 6, 9)) + 1),
                                   season != .) %>%
                            select(-season))) %>%
  unnest(col = c(train_data)) %>%
  # Nest by season and location
  nest(data = c(epiweek, ILI, year, week, region_hits, hits, order_week, month)) %>%
  # Create time series of ILI and numerical value of lambda
  mutate(data = map(data,
                   ~ mutate(.x,
                            ILI = ts(ILI, frequency = 52, start = c(2010, 40)))),
         lambda = 0,
         group = rep(1:length(cluster), length.out = nrow(.)))

# Split for parallel model fitting
fourier_fits_by_group <- fourier_model_fit_data %>%
  partition(group, cluster = cluster)

fourier_fits_by_group %>%
  cluster_library(c("tidyverse", "forecast", "lubridate", "FluSight", "MMWRweek"))

fourier_model_fits <- fourier_fits_by_group %>%
  # Fit model
  mutate(fit_1 = map(data,
                     ~ auto.arima(.x$ILI, xreg = fourier(.x$ILI, K = 1),
                                         seasonal = FALSE, lambda = 0)),
         fit_2 = map(data, 
                      ~ auto.arima(.x$ILI, xreg = fourier(.x$ILI, K = 2),
                                          seasonal = FALSE, lambda = 0)),
         fit_3 = map(data, 
                      ~ auto.arima(.x$ILI, xreg = fourier(.x$ILI, K = 3),
                                          seasonal = FALSE, lambda = 0)),
         fit_4 = map(data, 
                      ~ auto.arima(.x$ILI, xreg = fourier(.x$ILI, K = 4),
                                          seasonal = FALSE, lambda = 0)),
         fit_5 = map(data, 
                     ~ auto.arima(.x$ILI, xreg = fourier(.x$ILI, K = 5),
                                          seasonal = FALSE, lambda = 0)),
         fit_6 = map(data,
                      ~ auto.arima(.x$ILI, xreg = fourier(.x$ILI, K = 6),
                                          seasonal = FALSE, lambda = 0)),
         fit_7 = map(data,
                      ~ auto.arima(.x$ILI, xreg = fourier(.x$ILI, K = 7),
                                          seasonal = FALSE, lambda = 0)),
         fit_8 = map(data, 
                      ~ auto.arima(.x$ILI, xreg = fourier(.x$ILI, K = 8),
                                          seasonal = FALSE, lambda = 0)),
         fit_9 = map(data, 
                      ~ auto.arima(.x$ILI, xreg = fourier(.x$ILI, K = 9),
                                          seasonal = FALSE, lambda = 0)),
         fit_10 = map(data,
                       ~ auto.arima(.x$ILI, xreg = fourier(.x$ILI, K = 10),
                                           seasonal = FALSE, lambda = 0)),
         fit_11 = map(data, 
                       ~ auto.arima(.x$ILI, xreg = fourier(.x$ILI, K = 11),
                                           seasonal = FALSE, lambda = 0)),
         fit_12 = map(data, 
                       ~ auto.arima(.x$ILI, xreg = fourier(.x$ILI, K = 12),
                                           seasonal = FALSE, lambda = 0))) %>%
  collect() %>%
  as_tibble() %>%
  ungroup() %>%
  select(-data) %>%
  gather(key = "model", value = "fit", fit_1:fit_12)

saveRDS(fourier_model_fits, file = "Data/state_fourier_fits.Rds")

# Set up data for model fitting in parallel
fourier_model_data <- crossing(model = c("fit_1", "fit_2", "fit_3", "fit_4",
                                         "fit_5", "fit_6", "fit_7", "fit_8",
                                         "fit_9", "fit_10", "fit_11", "fit_12"),
                       season = c("2014/2015", "2015/2016", "2016/2017", 
                                  "2017/2018", "2018/2019"),
                       week = c(43:71),
                       location = unique(state_flu_data_merge$location)) %>%
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
  inner_join(fourier_model_fits, by = c("season", "location", "model")) %>%
  inner_join(ili_backfill_densities, by = c("location", 'week')) %>%
  # Set up grouping for parallel
  mutate(group = rep(1:length(cluster), length.out = nrow(.)))


# Set up party_df and load necessary libraries and functions
# load("Data/fourier_scores.Rdata")
fourier_start_time <- Sys.time()

message(paste('Start time:', fourier_start_time))

fourier_scores <- tibble()

for(this_season in unique(fourier_model_data$season)[3:5]) {
 
  fourier_by_group <- fourier_model_data %>%
    filter(season == this_season) %>%
    partition(group, cluster = cluster)
  
  fourier_by_group %>%
    cluster_library(c("tidyverse", "forecast", "lubridate", "FluSight", "MMWRweek")) %>%
    cluster_assign_value("state_flu_data_merge", state_flu_data_merge) %>%
    cluster_assign_value("ili_init_pub_list", ili_init_pub_list) %>%
    cluster_assign_value("fit_to_forecast", fit_to_forecast) %>%
    cluster_assign_value("sample_predictive_trajectories_arima",
                         sample_predictive_trajectories_arima)
  
  # Create forecasts for Fourier terms
  fourier_forecasts <- fourier_by_group %>%
    mutate(
      pred_data = pmap(
        list(season, week, location, epiweek),
           ~ filter(state_flu_data_merge, year <= as.numeric(substr(..1, 6, 9)),
                    season != paste0(substr(..1, 6, 9), "/",
                                     as.numeric(substr(..1, 6, 9)) + 1),
                    season != ..1 | order_week %in% 40:..2,
                    location == ..3) %>%
             select(epiweek, location, ILI, season, week, order_week) %>%
             left_join(select(ili_init_pub_list[[paste(..4)]],
                              ILI, epiweek, location),
                       by = c("epiweek", "location")) %>%
             mutate(ILI = ifelse(is.na(ILI.y), ILI.x, ILI.y)) %>%
             select(-ILI.x, -ILI.y)
      ),
      
      # Calculate max week for each season
      max_week = ifelse(season == "2014/2015", 53, 52),
      
      # Fit models,
      pred_results = pmap(
        list(fit, pred_data, location, season, model, max_week, backfill),
        ~ fit_to_forecast(object = ..1,
                          k = as.numeric(str_extract(..5, "[1-9][0-9]|[0-9]")),
                          pred_data = ..2,
                          location = ..3,
                          season = ..4,
                          max_week = ..6,
                          backfill = ..7,
                          npaths = 250)
      )
    ) %>%
    collect() %>%
    as_tibble() %>%
    ungroup() %>%
    select(season, model, location, week, pred_results)
  
 # Normalize probabilities and score forecasts
  fourier_scores <- fourier_forecasts %>%
    mutate(pred_results = map2(pred_results, location,
                               ~ mutate(.x, location = .y) %>%
                                 normalize_probs())) %>%
    select(season, model, week, pred_results) %>%
    unnest(col = c(pred_results)) %>%
    nest(data = c(target, bin_start_incl, value, type, unit, bin_end_notincl, 
                  forecast_week, location)) %>%
    left_join(state_nested_truth, by = "season") %>%
    mutate(scores = map2(data, truth,
                         ~ score_entry(.x, .y))) %>%
    select(season, model, scores) %>%
    unnest(col = c(scores)) 

  saveRDS(fourier_scores, file = paste0("Data/state_CV_Fourier_Scores_",
                                       substr(this_season, 1, 4), ".Rds"))
  
  rm(fourier_by_group, fourier_forecasts, fourier_scores)
  
  message(paste('Season', this_season, 'finished at', Sys.time()))
}

fourier_scores <- readRDS('Data/state_CV_Fourier_Scores_2014.Rds') %>%
  bind_rows(readRDS('Data/state_CV_Fourier_Scores_2015.Rds')) %>%
  bind_rows(readRDS('Data/state_CV_Fourier_Scores_2016.Rds')) %>%
  bind_rows(readRDS('Data/state_CV_Fourier_Scores_2017.Rds')) %>%
  bind_rows(readRDS('Data/state_CV_Fourier_Scores_2018.Rds'))

# Determine best K value for each region
best_k_cv <- fourier_scores %>%
  group_by(location, model) %>%
  summarize(avg_score = mean(score)) %>%
  filter(avg_score == max(avg_score)) %>%
  mutate(K = as.numeric(str_extract(model, "[1-9][0-9]|[0-9]"))) %>%
  ungroup() %>%
  select(location, K)

saveRDS(best_k_cv, file = "Data/state_CV_Fourier_terms.Rds")

##### ARIMA structure for error terms #####

best_k_cv <- readRDS("Data/state_CV_Fourier_terms.Rds")

# Create training fits
arima_model_fits <- readRDS("Data/state_arima_fits.Rds")

arima_model_fit_data <- crossing(season = c("2014/2015", "2015/2016", "2016/2017", 
                                            "2017/2018", "2018/2019"),
                                arima_1 = 0:3,
                                arima_2 = 0:1,
                                arima_3 = 0:3) %>%
  mutate(train_data = map(season,
                          ~ filter(state_flu_data_merge, year <= as.numeric(substr(., 6, 9)),
                                   season != paste0(substr(., 6, 9), "/",
                                                    as.numeric(substr(., 6, 9)) + 1),
                                   season != .) %>%
                            select(-season))) %>%
  unnest(col = c(train_data)) %>%
  # Nest by season and location
  nest(data = c(epiweek, ILI, year, week, region_hits, hits, order_week, month)) %>%
  # Create time series of ILI
  mutate(data = map(data,
                   ~ mutate(.x,
                            ILI = ts(ILI, frequency = 52, start = c(2010, 40))))) %>%
  # Join fourier values
  left_join(best_k_cv, by = "location") %>%
  # Set up for parallel
  mutate(group = rep(1:length(cluster), length.out = nrow(.)),
         lambda = 0) %>%
  select(-transform)

# Set up clusters
arima_fit_parallel <- arima_model_fit_data %>%
  partition(group, cluster = cluster)
arima_fit_parallel %>%
  cluster_library(c("tidyverse", "forecast", "lubridate", "FluSight", "MMWRweek"))

# Fit ARIMA models
arima_model_fits <- arima_fit_parallel %>%
  mutate(fit = pmap(
    list(data, arima_1, arima_2, arima_3, K, lambda),
    ~ try(Arima(..1$ILI, order = c(..2, ..3, ..4),
            xreg = fourier(..1$ILI, K = ..5),
            lambda = ..6), silent = TRUE)
    )) %>%
  select(-data) %>%
  collect() %>%
  as_tibble() %>%
  ungroup() %>%
  filter(!grepl("Error", fit))

saveRDS(arima_model_fits, file = "Data/state_arima_fits.Rds")

# Set up data for forecast creation in parallel
arima_model_data_setup <- crossing(season = c("2014/2015", "2015/2016", "2016/2017", 
                                              "2017/2018", "2018/2019"),
                             arima_1 = 0:3,
                             arima_2 = 0:1,
                             arima_3 = 0:3,
                             week = c(43:71),
                             location = unique(state_flu_data_merge$location)) %>%
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
  inner_join(select(arima_model_fits, -group), 
            by = c("location", "season", "arima_1", "arima_2", "arima_3")) %>%
  # Set up grouping for parallel
  mutate(group = rep(1:length(cluster), length.out = nrow(.)))

# Create and score forecasts in single-season parallel chunks
# Clear forecasts from memory after each season

start_time <- Sys.time()
for(this_season in unique(arima_model_data_setup$season)) {
  
  arima_model_data_parallel <- arima_model_data_setup %>%
    filter(season == this_season) %>%
    partition(group, cluster = cluster)
  
  arima_model_data_parallel %>%
    cluster_library(c("tidyverse", "forecast", "lubridate", "FluSight", "MMWRweek")) %>% 
    cluster_assign_value("state_flu_data_merge", state_flu_data_merge) %>%
    cluster_assign_value("ili_init_pub_list", ili_init_pub_list) %>%
    cluster_assign_value("fit_to_forecast", fit_to_forecast) %>%
    cluster_assign_value("sample_predictive_trajectories_arima", 
                         sample_predictive_trajectories_arima)
  

  # Create forecasts for different ARIMA structures
  arima_forecasts <- arima_model_data_parallel %>%
    mutate(pred_data = pmap(list(season, week, location, epiweek), 
                     ~ filter(state_flu_data_merge, year <= as.numeric(substr(..1, 6, 9)),
                              season != paste0(substr(..1, 6, 9), "/",
                                               as.numeric(substr(..1, 6, 9)) + 1),
                              season != ..1 | order_week %in% 40:..2,
                              location == ..3) %>%
                       left_join(select(ili_init_pub_list[[paste(..4)]], 
                                        ILI, epiweek, location),
                                 by = c("epiweek", "location")) %>%
                       mutate(ILI = ifelse(is.na(ILI.y), ILI.x, ILI.y)) %>%
                       select(-ILI.x, -ILI.y) %>%
                       mutate(ILI = ts(ILI, frequency = 52, start = c(2010, 40))))) %>%
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
                                       h = ..5 - 14 - nrow(..2[..2$season == ..4, ])),
                        pred_data = ..2,
                        location = ..3,
                        season = ..4,
                        max_week = ..5,
                        npaths = 250)
      )) %>%
    collect() %>%
    as_tibble() %>%
    ungroup() %>%
    select(season, arima_1:arima_3, location, week, pred_results)
 
  # Normalize probabilities and score forecasts 
  arima_scores <- arima_forecasts %>%
    mutate(pred_results = map2(pred_results, location,
                               ~ mutate(.x, location = .y) %>%
                                 normalize_probs())) %>%
    select(season, arima_1:arima_3, week, pred_results) %>%
    unnest(col = c(pred_results)) %>%
    nest(data = c(target, bin_start_incl, value, type, unit, bin_end_notincl, 
                  forecast_week, location)) %>%
    left_join(state_nested_truth, by = "season") %>%
    mutate(scores = map2(data, exp_truth,
                         ~ score_entry(.x, .y))) %>%
    select(season, arima_1:arima_3, scores) %>%
    unnest(col = c(scores)) 
  
  saveRDS(arima_scores, file = paste0("Data/state_CV_ARIMA_Scores_", 
                                      substr(this_season, 1, 4), ".Rds"))
  
  rm(arima_forecasts, arima_model_data_parallel)
  
  message(paste("Total elapsed time", lubridate::time_length(Sys.time() - start_time, unit = "min")))

}
arima_time <- lubridate::time_length(Sys.time() - start_time, unit = "min")

arima_scores_1415 <- readRDS("Data/state_CV_ARIMA_Scores_2014.Rds")
arima_scores_1516 <- readRDS("Data/state_CV_ARIMA_Scores_2015.Rds")
arima_scores_1617 <- readRDS("Data/state_CV_ARIMA_Scores_2016.Rds")
arima_scores_1718 <- readRDS("Data/state_CV_ARIMA_Scores_2017.Rds")
arima_scores_1819 <- readRDS("Data/state_CV_ARIMA_Scores_2018.Rds")

# Determine best ARIMA model for each region
best_arima_cv <-  bind_rows(arima_scores_1415, arima_scores_1516,
                            arima_scores_1617, arima_scores_1718,
                            arima_scores_1819) %>%
  group_by(location, arima_1, arima_2, arima_3) %>%
  summarize(avg_score = mean(score)) %>%
  group_by(location) %>%
  filter(avg_score == max(avg_score)) %>%
  ungroup() %>%
  select(location, arima_1:arima_3)

saveRDS(best_arima_cv, file = "Data/state_CV_ARIMA_terms.Rds")

##### Additional covariates #####
best_k_cv <- readRDS("data/state_CV_Fourier_terms.Rds")
best_arima_cv <- readRDS("Data/state_CV_ARIMA_terms.Rds")

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

covar_model_fits <- readRDS("Data/state_covar_fits.Rds")
covar_model_fits_data <- crossing(season = c("2014/2015", "2015/2016", "2016/2017", 
                                             "2017/2018", "2018/2019"),
                             model = c("ARIMA only", "National Gtrends", "FluVirus",
                                       "Backfill", "National Gtrends & Backfill",
                                       "FluVirus & Backfill",
                                       "National Gtrends & FluVirus",
                                       "National Gtrends, FluVirus, & Backfill",
                                       "Regional Gtrends", "Regional Gtrends & Backfill",
                                       "Regional Gtrends & FluVirus",
                                       "Regional Gtrends, FluVirus, & Backfill")) %>%
  mutate(train_data = map(season,
                          ~ filter(state_flu_data_merge, year <= as.numeric(substr(., 6, 9)),
                                   season != paste0(substr(., 6, 9), "/",
                                                    as.numeric(substr(., 6, 9)) + 1),
                                   season != .))) %>%
  unnest(col = c(train_data)) %>%
  # Nest by season and location
  nest(-season, -location, -model) %>%
  # Create time series of ILI
  mutate(data = map(data,
                   ~ mutate(.x,
                            ILI = ts(ILI, frequency = 52, start = c(2010, 40)))),
         lambda = 0,
         group = rep(1:length(cluster), length.out = nrow(.))) %>%
  # Merge best lambda, Fourier K value, and ARIMA structure by location
  left_join(best_k_cv, by = "location") %>%
  left_join(best_arima_cv, by = "location") %>%
# Remove virus models for states w/o virus data
filter(!(model %in% c("FluVirus", "FluVirus & Backfill", "National Gtrends & FluVirus",
                    "National Gtrends, FluVirus, & Backfill", "Regional Gtrends & FluVirus",
                    "Regional Gtrends, FluVirus, & Backfill") &
         location %in% c("Alaksa", "New Jersey", "New Hampshire", "Rhode Island", "Louisiana",
                         "Wyoming")),
       !(model %in% c("Backfill", "National Gtrends & Backfill",
                     "FluVirus & Backfill","National Gtrends, FluVirus, & Backfill",
                     "Regional Gtrends & Backfill", "Regional Gtrends, FluVirus, & Backfill") &
          location  %in% c("Oklahoma")),
       (location != "Louisiana" | season != "2017-2018"))

out_fits <- list()

for (this_location in unique(covar_model_fits_data$location)) {

  out_fits[[this_location]] <- filter(covar_model_fits_data, location == this_location) %>%
    # Fit models
mutate(fit = ifelse(
  model == "ARIMA only",
    pmap(
      list(data, arima_1, arima_2, arima_3, K, lambda),
      ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
              xreg = fourier(..1$ILI, K = ..5),
              lambda = ..6)
    ),
  ifelse(
    model == "National Gtrends",
      pmap(
        list(data, arima_1, arima_2, arima_3, K, lambda),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$hits) %>%
                  as.matrix(),
                lambda = ..6)
      ),
  ifelse(
    model == "Regional Gtrends",
      pmap(
        list(data, arima_1, arima_2, arima_3, K, lambda),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$region_hits)%>%
                  as.matrix(),
                lambda = ..6)
      ),
  ifelse(
    model == "FluVirus",
      pmap(
        list(data, arima_1, arima_2, arima_3, K, lambda),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$cum_h1per, ..1$cum_h3per)%>%
                  as.matrix(),
                lambda = ..6)
      ),
  ifelse(
    model == "Backfill",
      pmap(
        list(data, arima_1, arima_2, arima_3, K, lambda),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$backfill) %>%
                  as.matrix(),
                lambda = ..6)
      ),
  ifelse(
    model == "National Gtrends & Backfill",
     pmap(
        list(data, arima_1, arima_2, arima_3, K, lambda),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$hits, ..1$backfill)%>%
                  as.matrix(),
                lambda = ..6)
      ),
  ifelse(
    model == "Regional Gtrends & Backfill",
      pmap(
        list(data, arima_1, arima_2, arima_3, K, lambda),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$region_hits, ..1$backfill) %>%
                  as.matrix(),
                lambda = ..6)
      ),
  ifelse(
    model == "FluVirus & Backfill",
      pmap(
        list(data, arima_1, arima_2, arima_3, K, lambda),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$cum_h1per, ..1$cum_h3per,
                             ..1$backfill) %>%
                  as.matrix(),
                lambda = ..6)
      ),
  ifelse(
    model == "National Gtrends & FluVirus",
      pmap(
        list(data, arima_1, arima_2, arima_3, K, lambda),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$hits, ..1$cum_h1per,
                             ..1$cum_h3per) %>%
                  as.matrix(),
                lambda = ..6)
      ),
  ifelse(
    model == "Regional Gtrends & FluVirus",
      pmap(
        list(data, arima_1, arima_2, arima_3, K, lambda),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$region_hits, ..1$cum_h1per,
                             ..1$cum_h3per) %>%
                  as.matrix(),
                lambda = ..6)
      ),
  ifelse(
    model == "National Gtrends, FluVirus, & Backfill",
      pmap(
        list(data, arima_1, arima_2, arima_3, K, lambda),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$hits, ..1$cum_h1per,
                             ..1$cum_h3per,
                             ..1$backfill) %>%
                  as.matrix(),
                lambda = ..6)
      ),
  ifelse(
    model == "Regional Gtrends, FluVirus, & Backfill",
      pmap(
        list(data, arima_1, arima_2, arima_3, K, lambda),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$region_hits, ..1$cum_h1per,
                             ..1$cum_h3per,
                             ..1$backfill) %>%
                  as.matrix(),
                lambda = ..6)
      ),
    list()
  )))))))))))) # Close ridiculous number of ifelse statements
)
  message(paste(this_location, "finished"))
}


covar_model_fits <- bind_rows(out_fits)

saveRDS(covar_model_fits, file = "Data/state_covar_fits.Rds")

# Set up data for forecast creation in parallel
covar_model_data_setup <- crossing(season = c("2014/2015", "2015/2016", "2016/2017",
                                              "2017/2018", "2018/2019"),
                                   model = c("ARIMA only", "National Gtrends", "FluVirus", 
                                             "Backfill", "National Gtrends & Backfill",
                                             "FluVirus & Backfill",
                                             "National Gtrends & FluVirus",
                                             "National Gtrends, FluVirus, & Backfill",
                                             "Regional Gtrends", "Regional Gtrends & Backfill",
                                             "Regional Gtrends & FluVirus",
                                             "Regional Gtrends, FluVirus, & Backfill"),
                                   week = c(43:71),
                                   location = unique(state_flu_data_merge$location)) %>%
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
  mutate(group = rep(1:length(cluster), length.out = nrow(.)))

start_time <- Sys.time()

for (this_season in unique(covar_model_data_setup$season)) {
    
  # Set up dataset for forecasting in parallel
  covar_model_parallel <- covar_model_data_setup %>%
    filter(season == this_season) %>%
    partition(group, cluster = cluster)
  
  covar_model_parallel %>%
    cluster_library(c("tidyverse", "forecast", "lubridate", "FluSight", "MMWRweek")) %>%
    cluster_assign_value("state_flu_data_merge", state_flu_data_merge) %>%
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
                       ~ filter(state_flu_data_merge, year <= as.numeric(substr(..1, 6, 9)),
                                season != paste0(substr(..1, 6, 9), "/",
                                                 as.numeric(substr(..1, 6, 9)) + 1),
                                season != ..1 | order_week %in% 40:..2,
                                location == ..3) %>%
                         left_join(select(ili_init_pub_list[[paste(..4)]],
                                          ILI, epiweek, location),
                                   by = c("epiweek", "location")) %>%
                         mutate(ILI = ifelse(is.na(ILI.y), ILI.x, ILI.y)) %>%
                         select(-ILI.x, -ILI.y) %>%
                         mutate(ILI = ts(ILI, frequency = 52, start = c(2010, 40)))),
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
                       ..1$cum_h1per, ..1$cum_h3per)),
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
                       ..1$cum_h1per, ..1$cum_h3per,
                       ..1$backfill)),
        model == "National Gtrends & FluVirus" ~ 
          pmap(list(pred_data, K),
               ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2)),
                     ..1$hits, ..1$cum_h1per, 
                     ..1$cum_h3per)),
        model == "Regional Gtrends & FluVirus" ~ 
          pmap(list(pred_data, K),
               ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2)),
                       ..1$region_hits, ..1$cum_h1per, 
                       ..1$cum_h3per)),
        model == "National Gtrends, FluVirus, & Backfill" ~ 
          pmap(list(pred_data, K),
               ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2)),
                       ..1$hits, ..1$cum_h1per, 
                       ..1$cum_h3per, ..1$backfill)),
        model == "Regional Gtrends, FluVirus, & Backfill" ~ 
          pmap(list(pred_data, K),
               ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2)),
                       ..1$region_hits, ..1$cum_h1per, 
                       ..1$cum_h3per, ..1$backfill))
        ),
      xreg = map(xreg,
                 ~ as.matrix(.)),
      max_week = ifelse(season == "2014/2015", 53, 52),
      # Fit models
      pred_fit = pmap(list(pred_data, fit, xreg),
                      ~ Arima(..1$ILI, xreg = ..3, model = ..2)),
      # Create data frame of xreg terms for forecasting
      gtrend_forecast = pmap(
        list(pred_data, location, season, prev_season, week, max_week),
        ~ tibble(hits = c(state_flu_data_merge %>%
                            filter(location == ..2, season == ..3, 
                                   order_week == ..5 + 1) %>%
                            pull(hits),
                          state_flu_data_merge %>%
                            filter(location == ..2, season == ..4, 
                                   order_week > ..5 + 1, 
                                   order_week < ..6 + 26) %>%
                            pull(hits) *
                            mean(state_flu_data_merge %>%
                              filter(location == ..2, season == ..3,
                                     order_week <= ..5 + 1) %>%
                              pull(hits) /
                              state_flu_data_merge %>%
                                filter(location == ..2, season == ..4,
                                       order_week <= ..5 + 1) %>%
                                pull(hits))))
      ),
      reg_gtrend_forecast = pmap(
        list(pred_data, location, season, prev_season, week, max_week),
        ~ tibble(region_hits = c(state_flu_data_merge %>%
                            filter(location == ..2, season == ..3, 
                                   order_week == ..5 + 1) %>%
                            pull(region_hits),
                          state_flu_data_merge %>%
                            filter(location == ..2, season == ..4, 
                                   order_week > ..5 + 1, 
                                   order_week < ..6 + 26) %>%
                            pull(region_hits) *
                            mean(state_flu_data_merge %>%
                                   filter(location == ..2, season == ..3,
                                          order_week <= ..5 + 1) %>%
                                   pull(region_hits) /
                                   state_flu_data_merge %>%
                                   filter(location == ..2, season == ..4,
                                          order_week <= ..5 + 1) %>%
                                   pull(region_hits))))
      ),
      h1_per_forecast = pmap(
        list(pred_data, season, max_week),
        ~ rep(last(..1$cum_h1per), ..3 - 14 - 
                nrow(..1[..1$season == ..2, ]))
      ),
      h3_per_forecast = pmap(
        list(pred_data, season, max_week),
        ~ rep(last(..1$cum_h3per), ..3 - 14 - 
                nrow(..1[..1$season == ..2, ]))
      ),
      backfill = pmap(
        list(data, week, max_week, season),
        function(data, week, max_week, season) {
          out <- numeric()
          for(i in week:(max_week + 24)) {
            temp <- data[data$order_week == i, ]
            
            out <- c(out, predict(lm(backfill ~ year, 
                                     data = temp),
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
                    h1_per_forecast, h3_per_forecast),
               ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                             h = ..3 - 14 - 
                                               nrow(..1[..1$season == ..4, ]))),
                       data.frame(h1_per_samples = ..5,
                                  h3_per_samples = ..6)) %>%
                 rename(`..1$cum_h1per` = h1_per_samples,
                        `..1$cum_h3per` = h3_per_samples)),
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
                    h1_per_forecast, h3_per_forecast, backfill),
               ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                             h = ..3 - 14 - 
                                               nrow(..1[..1$season == ..4, ]))),
                       data.frame(h1_per_samples = ..5,
                                  h3_per_samples = ..6, 
                                  backfill = ..7)) %>%
                 rename(`..1$cum_h1per` = h1_per_samples,
                        `..1$cum_h3per` = h3_per_samples,
                        `..1$backfill` = backfill)),
        model == "National Gtrends & FluVirus" ~
          pmap(list(pred_data, K, max_week, season, gtrend_forecast,
                    h1_per_forecast, h3_per_forecast),
               ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                             h = ..3 - 14 - 
                                               nrow(..1[..1$season == ..4, ]))),
                       data.frame(hits = ..5, 
                                  h1_per_samples = ..6,
                                  h3_per_samples = ..7)) %>%
                 rename(`..1$hits` = hits,
                        `..1$cum_h1per` = h1_per_samples,
                        `..1$cum_h3per` = h3_per_samples)),
        model == "Regional Gtrends & FluVirus" ~
          pmap(list(pred_data, K, max_week, season, reg_gtrend_forecast,
                    h1_per_forecast, h3_per_forecast),
               ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                             h = ..3 - 14 - 
                                               nrow(..1[..1$season == ..4, ]))),
                       data.frame(region_hits = ..5, 
                                  h1_per_samples = ..6,
                                  h3_per_samples = ..7)) %>%
                 rename(`..1$region_hits` = region_hits,
                        `..1$cum_h1per` = h1_per_samples,
                        `..1$cum_h3per` = h3_per_samples)),
        model == "National Gtrends, FluVirus, & Backfill" ~
          pmap(list(pred_data, K, max_week, season, gtrend_forecast,
                    h1_per_forecast, h3_per_forecast, backfill),
               ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                             h = ..3 - 14 - 
                                               nrow(..1[..1$season == ..4, ]))),
                       data.frame(hits = ..5, 
                                  h1_per_samples = ..6,
                                  h3_per_samples = ..7, 
                                  backfill = ..8)) %>%
                 rename(`..1$hits` = hits,
                        `..1$cum_h1per` = h1_per_samples,
                        `..1$cum_h3per` = h3_per_samples,
                        `..1$backfill` = backfill)),
        model == "Regional Gtrends, FluVirus, & Backfill" ~
          pmap(list(pred_data, K, max_week, season, reg_gtrend_forecast,
                    h1_per_forecast, h3_per_forecast, backfill),
               ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                             h = ..3 - 14 - 
                                               nrow(..1[..1$season == ..4, ]))),
                       data.frame(region_hits = ..5, 
                                  h1_per_samples = ..6,
                                  h3_per_samples = ..7, 
                                  backfill = ..8)) %>%
                 rename(`..1$region_hits` = region_hits,
                        `..1$cum_h1per` = h1_per_samples,
                        `..1$cum_h3per` = h3_per_samples,
                        `..1$backfill` = backfill))
      ),
      forecast_xreg = map(forecast_xreg,
                          ~ as.matrix(.)),
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
    as_tibble()
  

  covar_scores <- covar_forecasts %>% ungroup() %>%
    select(season, model, location, week, pred_results) %>%
    mutate(pred_results = map2(pred_results, location,
                               ~ mutate(.x, location = .y) %>%
                                 normalize_probs())) %>%
    select(season, model, week, pred_results) %>%
    unnest(col = c(pred_results)) %>%
    nest(data = c(target, bin_start_incl, value, type, unit, bin_end_notincl, 
                  forecast_week, location)) %>%
    left_join(state_nested_truth, by = "season") %>%
    mutate(scores = map2(data, truth,
                         ~ score_entry(.x, .y))) %>%
    select(season, model, scores) %>%
    unnest(co = c(scores)) 
  
  saveRDS(covar_scores, file = paste0("Data/state_CV_covar_Scores_", 
                                      substr(this_season, 1, 4), ".Rds"))
  
  message(paste(this_season, "season complete. Total elapsed time:", 
                Sys.time() - start_time))
  
  rm(covar_forecasts, covar_model_parallel, cluster)
  
}

covar_scores_1415 <- readRDS("Data/state_CV_covar_Scores_2014.Rds")
covar_scores_1516 <- readRDS("Data/state_CV_covar_Scores_2015.Rds")
covar_scores_1617 <- readRDS("Data/state_CV_covar_Scores_2016.Rds")
covar_scores_1718 <- readRDS("Data/state_CV_covar_Scores_2017.Rds")
covar_scores_1819 <- readRDS("Data/state_CV_covar_Scores_2018.Rds")

# Determine best covariates for each region
best_covar_cv <- bind_rows(covar_scores_1415, covar_scores_1516,
                           covar_scores_1617, covar_scores_1718,
                           covar_scores_1819) %>%
  group_by(location, model) %>%
  summarize(avg_score = mean(score)) %>%
  group_by(location) %>%
  filter(avg_score == max(avg_score)) %>%
  ungroup() %>%
  select(location, model)

# Manually select model for each region from results
# If multiple models similar (~ 0.01 diff), take model with most information
# best_covar_cv <- tibble(location = c("HHS Region 1", "HHS Region 2", "HHS Region 3",
#                                      "HHS Region 4", "HHS Region 5", "HHS Region 6",
#                                      "HHS Region 7", "HHS Region 8", "HHS Region 9",
#                                      "HHS Region 10", "US National"),
#                         model = c("Regional Gtrends", "FluVirus & Backfill",
#                                   "Regional Gtrends", "National Gtrends & Backfill",
#                                   "FluVirus", "ARIMA only",
#                                   "Regional Gtrends & Backfill",
#                                   "Regional Gtrends, FluVirus, & Backfill",
#                                   "Regional Gtrends & Backfill",
#                                   "Regional Gtrends",
#                                   "National Gtrends & FluVirus"))

saveRDS(best_covar_cv, file = "Data/state_CV_covar_terms.Rds")


##### Create final fits and build forecasts #####
best_k_cv <- readRDS("Data/state_CV_Fourier_terms.Rds")
best_arima_cv <- readRDS("Data/state_CV_ARIMA_terms.Rds")
best_covar_cv <- readRDS("Data/state_CV_covar_terms.Rds")

# Set up model fits
final_fits <- readRDS("Data/state_Final_CV_fits.Rds")
final_fits_data <- tibble(season = c("2014/2015", "2015/2016",
                                     "2016/2017", "2017/2018",
                                     "2018/2019")) %>%
  mutate(train_data = map(season,
                          ~ filter(state_flu_data_merge, year <= as.numeric(substr(., 6, 9)),
                                   season != paste0(substr(., 6, 9), "/",
                                                    as.numeric(substr(., 6, 9)) + 1),
                                   season != .) %>%
                            select(-season))) %>%
  unnest(col = c(train_data)) %>%
  # Nest by season and location
  nest(data = c(epiweek, ILI, year, week, region_hits, hits, orig_ILI, final_ILI, 
                backfill, order_week, month)) %>%
  # Create time series of ILI
  mutate(data = map(data,
                    ~ mutate(.x,
                             ILI = ts(ILI, frequency = 52, start = c(2010, 40)))),
         lambda = 0) %>%
  # Merge best lambda, Fourier K value, and ARIMA structure by location
  left_join(best_k_cv, by = "location") %>%
  left_join(best_arima_cv, by = "location") %>%
  left_join(best_covar_cv, by = "location") %>%
  # Remove virus models for states w/o virus data
  filter(!(model %in% c("FluVirus", "FluVirus & Backfill", "National Gtrends & FluVirus",
                        "National Gtrends, FluVirus, & Backfill", "Regional Gtrends & FluVirus",
                        "Regional Gtrends, FluVirus, & Backfill") &
             location %in% c("Alaksa", "New Jersey", "New Hampshire", "Rhode Island", "Louisiana",
                             "Wyoming")),
         !(model %in% c("Backfill", "National Gtrends & Backfill",
                        "FluVirus & Backfill","National Gtrends, FluVirus, & Backfill",
                        "Regional Gtrends & Backfill", "Regional Gtrends, FluVirus, & Backfill") &
             location  %in% c("Oklahoma")),
         (location != "Louisiana" | season != "2017-2018")) #%>%
  # Fit models

final_fit_list <- list()

for (this_location in unique(final_fits_data$location)) {
  final_fit_list[[this_location]] <- final_fits_data %>%
    filter(location == this_location) %>%
    mutate(fit = ifelse(
      model == "ARIMA only",
      pmap(
        list(data, arima_1, arima_2, arima_3, K, lambda),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = fourier(..1$ILI, K = ..5),
                lambda = ..6)
      ),
    ifelse(
      model == "National Gtrends",
      pmap(
        list(data, arima_1, arima_2, arima_3, K, lambda),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$hits) %>%
                  as.matrix(),
                lambda = ..6)
    ),
    ifelse(
      model == "Regional Gtrends",
      pmap(
        list(data, arima_1, arima_2, arima_3, K, lambda),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$region_hits)%>%
                  as.matrix(),
                lambda = ..6)
    ),
    ifelse(
      model == "FluVirus",
      pmap(
        list(data, arima_1, arima_2, arima_3, K, lambda),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$cum_h1per, ..1$cum_h3per)%>%
                  as.matrix(),
                lambda = ..6)
    ),
    ifelse(
      model == "Backfill",
      pmap(
        list(data, arima_1, arima_2, arima_3, K, lambda),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$backfill) %>%
                  as.matrix(),
                lambda = ..6)
    ),
    ifelse(
      model == "National Gtrends & Backfill",
      pmap(
        list(data, arima_1, arima_2, arima_3, K, lambda),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$hits, ..1$backfill)%>%
                  as.matrix(),
                lambda = ..6)
    ),
    ifelse(
      model == "Regional Gtrends & Backfill",
      pmap(
        list(data, arima_1, arima_2, arima_3, K, lambda),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$region_hits, ..1$backfill) %>%
                  as.matrix(),
                lambda = ..6)
    ),
    ifelse(
      model == "FluVirus & Backfill",
      pmap(
        list(data, arima_1, arima_2, arima_3, K, lambda),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$cum_h1per, ..1$cum_h3per,
                             ..1$backfill) %>%
                  as.matrix(),
                lambda = ..6)
    ),
    ifelse(
      model == "National Gtrends & FluVirus",
      pmap(
        list(data, arima_1, arima_2, arima_3, K, lambda),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$hits, ..1$cum_h1per,
                             ..1$cum_h3per) %>%
                  as.matrix(),
                lambda = ..6)
    ),
    ifelse(
      model == "Regional Gtrends & FluVirus",
      pmap(
        list(data, arima_1, arima_2, arima_3, K, lambda),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$region_hits, ..1$cum_h1per,
                             ..1$cum_h3per) %>%
                  as.matrix(),
                lambda = ..6)
    ),
    ifelse(
      model == "National Gtrends, FluVirus, & Backfill",
      pmap(
        list(data, arima_1, arima_2, arima_3, K, lambda),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$hits, ..1$cum_h1per,
                             ..1$cum_h3per,
                             ..1$backfill) %>%
                  as.matrix(),
                lambda = ..6)
    ),
    ifelse(
      model == "Regional Gtrends, FluVirus, & Backfill",
      pmap(
        list(data, arima_1, arima_2, arima_3, K, lambda),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$region_hits, ..1$cum_h1per,
                             ..1$cum_h3per,
                             ..1$backfill) %>%
                  as.matrix(),
                lambda = ..6)
    ),
    list()
    )))))))))))) # Close ridiculous number of ifelse statements
  )

  message(paste(this_location, "fits complete"))
}

final_fits <- bind_rows(final_fit_list)

saveRDS(final_fits, file = "Data/state_Final_CV_fits.Rds")

final_forecast_data_setup <- crossing(season = c("2014/2015", "2015/2016",
                                             "2016/2017", "2017/2018",
                                             "2018/2019"),
                                      week = c(40:73),
                                      location = unique(state_flu_data_merge$location)) %>%
  filter((week < 73 | season == "2014/2015") & 
         (season != "2018/2019" | week <= 70)) %>%
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
    prev_season = paste0(as.numeric(substr(season, 1, 4)) - 1, "/",
                    substr(season, 1, 4))
  ) %>%
  inner_join(final_fits, 
             by = c("location", "season")) 

# Create and save forecast files

for(this_season in unique(final_forecast_data_setup$season)) {
  
  # Create folder if needed
  season_path <- paste0("State Forecasts/", substr(this_season, 1, 4), "-",
                        substr(this_season, 6, 9), "/Dynamic Harmonic Model/")

  dir.create(season_path, showWarnings = FALSE)
  
  for(this_week in unique(final_forecast_data_setup$week[final_forecast_data_setup$season == this_season])) {

    temp <- final_forecast_data_setup %>%
      filter(season == this_season, week == this_week) %>%
      mutate(
        pred_data = pmap(list(season, week, location, epiweek), 
                         ~ filter(state_flu_data_merge, year <= as.numeric(substr(..1, 6, 9)),
                                  season != paste0(substr(..1, 6, 9), "/",
                                                   as.numeric(substr(..1, 6, 9)) + 1),
                                  season != ..1 | order_week %in% 40:..2,
                                  location == ..3) %>%
                           left_join(select(ili_init_pub_list[[paste(..4)]],
                                            ILI, epiweek, location),
                                     by = c("epiweek", "location")) %>%
                           mutate(ILI = ifelse(is.na(ILI.y), ILI.x, ILI.y)) %>%
                           select(-ILI.x, -ILI.y) %>%
                           mutate(ILI = ts(ILI, frequency = 52, start = c(2010, 40)))),
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
                         ..1$cum_h1per, ..1$cum_h3per)),
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
                         ..1$cum_h1per, ..1$cum_h3per,
                         ..1$backfill)),
          model == "National Gtrends & FluVirus" ~ 
            pmap(list(pred_data, K),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2)),
                         ..1$hits, ..1$cum_h1per, 
                         ..1$cum_h3per)),
          model == "Regional Gtrends & FluVirus" ~ 
            pmap(list(pred_data, K),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2)),
                         ..1$region_hits, ..1$cum_h1per, 
                         ..1$cum_h3per)),
          model == "National Gtrends, FluVirus, & Backfill" ~ 
            pmap(list(pred_data, K),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2)),
                         ..1$hits, ..1$cum_h1per, 
                         ..1$cum_h3per, ..1$backfill)),
          model == "Regional Gtrends, FluVirus, & Backfill" ~ 
            pmap(list(pred_data, K),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2)),
                         ..1$region_hits, ..1$cum_h1per, 
                         ..1$cum_h3per, ..1$backfill))
        ),
        xreg = map(xreg,
                   ~ as.matrix(.)),
        max_week = ifelse(season == "2014/2015", 53, 52),
        # Fit models
        pred_fit = pmap(list(pred_data, fit, xreg),
                        ~ Arima(..1$ILI, xreg = ..3, model = ..2)),
        # Create data frame of xreg terms for forecasting
        gtrend_forecast = pmap(
          list(pred_data, location, season, prev_season, week, max_week),
          ~ tibble(hits = c(state_flu_data_merge %>%
                              filter(location == ..2, season == ..3, 
                                     order_week == ..5 + 1) %>%
                              pull(hits),
                            state_flu_data_merge %>%
                              filter(location == ..2, season == ..4, 
                                     order_week > ..5 + 1, 
                                     order_week < ..6 + 26) %>%
                              pull(hits) *
                              mean(state_flu_data_merge %>%
                                     filter(location == ..2, season == ..3,
                                            order_week <= ..5 + 1) %>%
                                     pull(hits) /
                                     state_flu_data_merge %>%
                                     filter(location == ..2, season == ..4,
                                            order_week <= ..5 + 1) %>%
                                     pull(hits))))
        ),
        reg_gtrend_forecast = pmap(
          list(pred_data, location, season, prev_season, week, max_week),
          ~ tibble(region_hits = c(state_flu_data_merge %>%
                                     filter(location == ..2, season == ..3, 
                                            order_week == ..5 + 1) %>%
                                     pull(region_hits),
                                   state_flu_data_merge %>%
                                     filter(location == ..2, season == ..4, 
                                            order_week > ..5 + 1, 
                                            order_week < ..6 + 26) %>%
                                     pull(region_hits) *
                                     mean(state_flu_data_merge %>%
                                            filter(location == ..2, season == ..3,
                                                   order_week <= ..5 + 1) %>%
                                            pull(region_hits) /
                                            state_flu_data_merge %>%
                                            filter(location == ..2, season == ..4,
                                                   order_week <= ..5 + 1) %>%
                                            pull(region_hits))))
        ),
        h1_per_forecast = pmap(
          list(pred_data, season, max_week),
          ~ rep(last(..1$cum_h1per), ..3 - 14 - 
                  nrow(..1[..1$season == ..2, ]))
        ),
        h3_per_forecast = pmap(
          list(pred_data, season, max_week),
          ~ rep(last(..1$cum_h3per), ..3 - 14 - 
                  nrow(..1[..1$season == ..2, ]))
        ),
        backfill = pmap(
          list(data, week, max_week, season),
          function(data, week, max_week, season) {
            out <- numeric()
            for(i in week:(max_week + 24)) {
              temp <- data[data$order_week == i, ]
              
              out <- c(out, predict(lm(backfill ~ year, 
                                       data = temp),
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
                      h1_per_forecast, h3_per_forecast),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                               h = ..3 - 14 - 
                                                 nrow(..1[..1$season == ..4, ]))),
                         data.frame(h1_per_samples = ..5,
                                    h3_per_samples = ..6)) %>%
                   rename(`..1$cum_h1per` = h1_per_samples,
                          `..1$cum_h3per` = h3_per_samples)),
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
                      h1_per_forecast, h3_per_forecast, backfill),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                               h = ..3 - 14 - 
                                                 nrow(..1[..1$season == ..4, ]))),
                         data.frame(h1_per_samples = ..5,
                                    h3_per_samples = ..6,
                                    backfill = ..7)) %>%
                   rename(`..1$cum_h1per` = h1_per_samples,
                          `..1$cum_h3per` = h3_per_samples,
                          `..1$backfill` = backfill)),
          model == "National Gtrends & FluVirus" ~
            pmap(list(pred_data, K, max_week, season, gtrend_forecast,
                      h1_per_forecast, h3_per_forecast),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                               h = ..3 - 14 - 
                                                 nrow(..1[..1$season == ..4, ]))),
                         data.frame(hits = ..5, 
                                    h1_per_samples = ..6,
                                    h3_per_samples = ..7)) %>%
                   rename(`..1$hits` = hits,
                          `..1$cum_h1per` = h1_per_samples,
                          `..1$cum_h3per` = h3_per_samples)),
          model == "Regional Gtrends & FluVirus" ~
            pmap(list(pred_data, K, max_week, season, reg_gtrend_forecast,
                      h1_per_forecast, h3_per_forecast),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                               h = ..3 - 14 - 
                                                 nrow(..1[..1$season == ..4, ]))),
                         data.frame(region_hits = ..5, 
                                    h1_per_samples = ..6,
                                    h3_per_samples = ..7)) %>%
                   rename(`..1$region_hits` = region_hits,
                          `..1$cum_h1per` = h1_per_samples,
                          `..1$cum_h3per` = h3_per_samples)),
          model == "National Gtrends, FluVirus, & Backfill" ~
            pmap(list(pred_data, K, max_week, season, gtrend_forecast,
                      h1_per_forecast, h3_per_forecast, backfill),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                               h = ..3 - 14 - 
                                                 nrow(..1[..1$season == ..4, ]))),
                         data.frame(hits = ..5, 
                                    h1_per_samples = ..6,
                                    h3_per_samples = ..7,
                                    backfill = ..8)) %>%
                   rename(`..1$hits` = hits,
                          `..1$cum_h1per` = h1_per_samples,
                          `..1$cum_h3per` = h3_per_samples,
                          `..1$backfill` = backfill)),
          model == "Regional Gtrends, FluVirus, & Backfill" ~
            pmap(list(pred_data, K, max_week, season, reg_gtrend_forecast,
                      h1_per_forecast, h3_per_forecast, backfill),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                               h = ..3 - 14 - 
                                                 nrow(..1[..1$season == ..4, ]))),
                         data.frame(region_hits = ..5, 
                                    h1_per_samples = ..6,
                                    h3_per_samples = ..7,
                                    backfill = ..8)) %>%
                   rename(`..1$region_hits` = region_hits,
                          `..1$cum_h1per` = h1_per_samples,
                          `..1$cum_h3per` = h3_per_samples,
                          `..1$backfill` = backfill))
        ),
        forecast_xreg = map(forecast_xreg,
                            ~ as.matrix(.)),
        pred_results = pmap(
          list(pred_fit, pred_data, forecast_xreg, location, season, max_week),
          ~ fit_to_forecast(object = ..1,
                            xreg = ..3,
                            pred_data = ..2,
                            location = ..4,
                            season = ..5,
                            max_week = ..6,
                            npaths = 500,
                            challenge = "state_ili")
        )
      ) %>%
      select(location, week, max_week, pred_results)
    
    # Save EW and year as objects
    EW <- case_when(temp$week[1] > temp$max_week[1] ~ 
                      str_pad(temp$week[1] - temp$max_week[1], 2, side = "left", 0),
                    TRUE ~ str_pad(temp$week[1], 2, "left", 0))
    
    this_year <- case_when(
      as.numeric(EW) >= 40 ~ substr(this_season, 1, 4),
      as.numeric(EW) < 40 ~ substr(this_season, 6, 9)
    )
    
    # Write CSVs to internal ensemble folders and FSN folder
    temp %>%
      unnest(col = c(pred_results)) %>%
      normalize_probs() %>%
      select(-week, -max_week, -forecast_week) %>%
      bind_rows(generate_point_forecasts(.)) %>%
      select(location, target, type, unit, bin_start_incl, bin_end_notincl, value) %>%
      write_csv(path = paste0(season_path, "/EW", EW, ".csv"))
    
  }
  
}

