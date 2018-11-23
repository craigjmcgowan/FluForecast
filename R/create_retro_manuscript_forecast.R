# Create weekly prospective forecasts for 2018-2019 season

library(tidyverse)
library(forecast)
library(lubridate)
library(FluSight)
library(MMWRweek)

source("R/utils.R")

# Load data
load("Data/ili.Rdata")
load("Data/virologic.Rdata")
load("Data/Gtrends.Rdata")
load("Data/state_gtrend.Rdata")

set.seed(4321) # For reproducibility of backfill sample

flu_data_merge <- select(ili_current, epiweek, ILI, year, week, season, location) %>%
  inner_join(select(virologic_combined, location, season, year, week, cum_h1per,
                   cum_h3per),
            by = c("location", "season", "year", "week")) %>%
  # Join backfill and replace missing backfill data with random draw from 
  # same season/week observed values
  full_join(select(ili_backfill, -orig_ILI, -final_ILI),
            by = c("location", "season", "year", "week")) %>%
  filter(!season %in% c("2003/2004", "2018/2019"), year >= 2004) %>%
  # Replace missing backfill data with Gaussian random draw from same season/week observed values
  mutate(
      sim_backfill = map2_dbl(location, week,
                              ~ rnorm(1,
                                 ili_backfill_avg$avg_backfill[
                                   ili_backfill_avg$location == .x &
                                     ili_backfill_avg$week == .y],
                                 ili_backfill_avg$sd_backfill[
                                   ili_backfill_avg$location == .x &
                                     ili_backfill_avg$week == .y])),
      backfill = case_when(
        !is.na(backfill) ~ backfill,
        TRUE ~ sim_backfill
        )
      ) %>%
  select(-sim_backfill) %>%
  # Add Google Trends data by location
  left_join(bind_rows(mutate(gtrend_US_flu_merge, location = "US National"),
                       mutate(gtrend_MA_flu_merge, location = "HHS Region 1"),
                       mutate(gtrend_NY_flu_merge, location = "HHS Region 2"),
                       mutate(gtrend_PA_flu_merge, location = "HHS Region 3"),
                       mutate(gtrend_FL_flu_merge, location = "HHS Region 4"),
                       mutate(gtrend_IL_flu_merge, location = "HHS Region 5"),
                       mutate(gtrend_TX_flu_merge, location = "HHS Region 6"),
                       mutate(gtrend_MO_flu_merge, location = "HHS Region 7"),
                       mutate(gtrend_CO_flu_merge, location = "HHS Region 8"),
                       mutate(gtrend_CA_flu_merge, location = "HHS Region 9"),
                       mutate(gtrend_WA_flu_merge, location = "HHS Region 10")) %>%
               rename(region_hits = hits) %>%
               select(-date),
             by = c("location", "season", "year", "week")) %>%
  # Add national Google Trends data
  left_join(gtrend_US_flu_merge,
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

nested_truth <- ili_current %>%
  filter(year >= 2010, !season %in% c("2009/2010", "2018/2019")) %>%
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
  

##### 2014-2015 #####
load("Data/CV_Transform_terms_1415.Rdata")
load("data/CV_Fourier_terms_1415.Rdata")
load("Data/CV_ARIMA_terms_1415.Rdata")
load("Data/CV_covar_terms_1415.Rdata")

fits <- filter(flu_data_merge, year <= 2014,
               season != "2014/2015") %>%
  # Nest by location
  nest(-location) %>%
  # Create time series of ILI
  mutate(data = map(data,
                    ~ mutate(.x,
                             ILI = ts(ILI, frequency = 52, start = c(2006, 40))))) %>%
  # Merge best lambda, Fourier K value, and ARIMA structure by location
  left_join(rename(best_transform_cv, transform = model), by = "location") %>%
  mutate(lambda = map_dbl(data,
                             ~ BoxCox.lambda(.$ILI)),
         lambda = case_when(transform == "no_trans" ~ NA_real_,
                            transform == "log" ~ 0,
                            transform == "Box_Cox" ~ lambda)) %>%
  left_join(best_k_cv, by = "location") %>%
  left_join(best_arima_cv, by = "location") %>%
  left_join(best_covar_cv, by = "location") %>%
  # Fit models
  mutate(fit = case_when(
    model == "ARIMA only" ~
      pmap(
        list(data, arima_1, arima_2, arima_3, K, lambda),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = fourier(..1$ILI, K = ..5),
                lambda = ..6)
      ),
    model == "National Gtrends" ~
      pmap(
        list(data, arima_1, arima_2, arima_3, K, lambda),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$hits),
                lambda = ..6)
      ),
    model == "Regional Gtrends" ~
      pmap(
        list(data, arima_1, arima_2, arima_3, K, lambda),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$region_hits),
                lambda = ..6)
      ),
    model == "FluVirus" ~
      pmap(
        list(data, arima_1, arima_2, arima_3, K, lambda),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$cum_h1per, ..1$cum_h3per),
                lambda = ..6)
      ),
    model == "Backfill" ~
      pmap(
        list(data, arima_1, arima_2, arima_3, K, lambda),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$backfill),
                lambda = ..6)
      ),
    model == "FluVirus & Backfill" ~
      pmap(
        list(data, arima_1, arima_2, arima_3, K, lambda),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$cum_h1per, ..1$cum_h3per,
                             ..1$backfill),
                lambda = ..6)
      ),
    model == "National Gtrends & FluVirus" ~
      pmap(
        list(data, arima_1, arima_2, arima_3, K, lambda),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$hits, ..1$cum_h1per,
                             ..1$cum_h3per),
                lambda = ..6)
      ),
    model == "Regional Gtrends & FluVirus" ~
      pmap(
        list(data, arima_1, arima_2, arima_3, K, lambda),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$region_hits, ..1$cum_h1per,
                             ..1$cum_h3per),
                lambda = ..6)
      ),
    model == "National Gtrends & Backfill" ~
      pmap(
        list(data, arima_1, arima_2, arima_3, K, lambda),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$hits, ..1$backfill),
                lambda = ..6)
      ),
    model == "Regional Gtrends & Backfill" ~
      pmap(
        list(data, arima_1, arima_2, arima_3, K, lambda),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$region_hits, ..1$backfill),
                lambda = ..6)
      ),
    model == "National Gtrends, FluVirus, & Backfill" ~
      pmap(
        list(data, arima_1, arima_2, arima_3, K, lambda),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$hits, ..1$cum_h1per,
                             ..1$cum_h3per, ..1$backfill),
                lambda = ..6)
      ),
    model == "Regional Gtrends, FluVirus, & Backfill" ~
      pmap(
        list(data, arima_1, arima_2, arima_3, K, lambda),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$region_hits, ..1$cum_h1per,
                             ..1$cum_h3per, ..1$backfill),
                lambda = ..6)
      )
  )) %>%
  select(-data)


forecasts_1415 <- crossing(season = "2014/2015",
                           week = 40:71,
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
  inner_join(fits, by = "location") %>%
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
      ~ rep(last(..1$cum_h1per), ..3 - 14 - 
              nrow(..1[..1$season == ..2, ]))
    ),
    h3_per_forecast = pmap(
      list(pred_data, season, max_week),
      ~ rep(last(..1$cum_h3per), ..3 - 14 - 
              nrow(..1[..1$season == ..2, ]))
    ),
    backfill = pmap(
      list(pred_data, week, max_week, season),
      function(data, week, max_week, season) {
        out <- numeric()
        for(i in week:(max_week + 24)) {
          
          temp <- data[data$order_week == i &
                           !data$season %in% c("2004/2005", "2005/2006", "2006/2007",
                                               "2008/2009", "2014/2015"), ]
          
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
    pred_results = pmap(
      list(pred_fit, pred_data, forecast_xreg, location, season, max_week),
      ~ fit_to_forecast(object = ..1,
                        xreg = ..3,
                        pred_data = ..2,
                        location = ..4,
                        season = ..5,
                        max_week = ..6,
                        npaths = 1000)
    )
  ) %>%
  select(location, week, pred_results) 

prospective_scores_1415 <- forecasts_1415 %>% 
  select(location, week, pred_results) %>%
  mutate(pred_results = map2(pred_results, location,
                             ~ mutate(.x, location = .y) %>%
                               normalize_probs()),
         season = "2014/2015") %>%
  select(-location) %>%
  unnest() %>%
  nest(-week, -season) %>%
  left_join(nested_truth, by = "season") %>%
  mutate(scores = map2(data, exp_truth,
                       ~ score_entry(.x, .y)),
         eval_scores = pmap(list(scores, eval_period, season),
                            ~ create_eval_scores(..1, ..2, ..3))) %>%
  select(season, eval_scores) %>%
  unnest() 

save(prospective_scores_1415, file = "Data/prospective_scores_1415.Rdata")
