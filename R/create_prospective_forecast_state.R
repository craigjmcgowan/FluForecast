# Create weekly prospective forecasts for 2018-2019 season

library(tidyverse)
library(FluSight)
library(MMWRweek)
library(cdcfluview)
library(forecast)
library(gtrendsR)
library(lubridate)
library(zoo)

source("R/utils.R")
source("R/create_subtype_forecast.R")

##### Set week that forecasts are being based on #####
EW <- 06
epiweek <- 201906

##### Update data #####
source("R/read_data_state.R")
source("R/gtrends_pull_data.R")

set.seed(101085) # For reproducibility of backfill sample

state_flu_data_merge <- select(state_current, epiweek, ILI, year, week, season, location) %>%
  inner_join(select(state_virologic, location, season, year, week, cum_h1per,
                   cum_h3per, cum_bper),
            by = c("location", "season", "year", "week")) %>%
  # Add Google Trends data by location
  right_join(bind_rows(mutate(gtrend_PA_flu_merge, location = "Pennsylvania"),
                       mutate(gtrend_DE_flu_merge, location = "Delaware"),
                       mutate(gtrend_MD_flu_merge, location = "Maryland"),
                       mutate(gtrend_VA_flu_merge, location = "Virginia"),
                       mutate(gtrend_WV_flu_merge, location = "West Virginia")) %>%
               filter(year > 2010 | (year == 2010 & week >= 40)) %>%
               rename(region_hits = hits) %>%
               select(-date),
            by = c("location", "season", "year", "week")) %>%
  right_join(filter(gtrend_US_flu_merge, year > 2010 | (year == 2010 & week >= 40)),
            by = c("season", "year", "week")) %>%
  full_join(state_backfill,
            by = c("location", "season", "year", "week")) %>%
  # Remove week 33 in 2014 so all seasons have 52 weeks - minimal activity
  filter(!(year == 2014 & week == 33)) %>%
  mutate(
    order_week = case_when(
      week < 40 & season == "2014/2015" ~ week + 53,
      week < 40 ~ week + 52,
      TRUE ~ week
    ),
    month = month(MMWRweek2Date(year, week))
  ) %>%
  # Replace missing backfill data with Gaussian random draw from same season/week observed values
  mutate(
    sim_backfill = map2_dbl(location, month,
                            ~ rnorm(1,
                                    state_backfill_month_avg$avg_backfill[
                                      state_backfill_month_avg$location == .x &
                                        state_backfill_month_avg$month == .y],
                                    state_backfill_month_avg$sd_backfill[
                                      state_backfill_month_avg$location == .x &
                                        state_backfill_month_avg$month == .y])),
      backfill = case_when(
        !is.na(backfill) ~ backfill,
        TRUE ~ sim_backfill
        )
      ) %>%
  select(-sim_backfill)
  

##### Kudu #####
vir_ssn_per <- state_virologic %>%
  group_by(season, location) %>%
  summarize(h1per = last(cum_h1per),
            h3per = last(cum_h3per)) 

# Create forecasts 
kudu_ili <- state_current %>%
  filter(year <= 2018, season != "2018/2019")

kudu_virologic <- state_virologic %>%
  filter(season == "2018/2019") %>%
  # Create variable to keep weeks in order
  mutate(order_week = week_inorder(week, season)) %>%
  select(year,location, week, h1per = h1per_of_a, 
         h3per = h3per_of_a, order_week) %>%
  select(location, order_week, h1per, h3per)

# Create directory to store forecasts
dir.create("State Forecasts/2018-2019/Subtype Historical Average/",
           showWarnings = FALSE)

# Create target densities and functions
subtype_densities_1819 <- create_subtype_densities(
  kudu_ili, vir_ssn_per, challenge = "state_ili"
)

subtype_functions_1819 <- modify_depth(
  subtype_densities_1819, 3,
  function(dens) approxfun(dens$x, dens$y, rule = 2)
)

order_week <- ifelse(EW < 40, EW + 52, EW)

kudu_pred <- create_subtype_forecast(
    functions = subtype_functions_1819,
    virologic = kudu_virologic,
    pub_week = order_week,
    season = "2018/2019",
    challenge = "state_ili"
  ) %>%
  select(location, target, type, unit, bin_start_incl, bin_end_notincl, value)
  
EW_paste <- str_pad(EW, 2, pad = "0")
  
write_csv(kudu_pred,
          path = paste0("State Forecasts/2018-2019/Subtype Historical Average/EW",
                        EW_paste, ".csv"))

##### Springbok #####
load("Data/state_CV_Transform_terms.Rdata")
load("data/state_CV_Fourier_terms.Rdata")
load("Data/state_CV_ARIMA_terms.Rdata")
load("Data/state_CV_covar_terms.Rdata")

fits <- filter(state_flu_data_merge, year <= 2018,
               season != "2018/2019") %>%
  # Nest by location
  nest(-location) %>%
  # Create time series of ILI
  mutate(data = map(data,
                    ~ mutate(.x,
                             ILI = ts(ILI, frequency = 52, start = c(2006, 40))))) %>%
  # Merge best lambda, Fourier K value, and ARIMA structure by location
  left_join(best_transform_cv, by = "location") %>%
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
                             ..1$hits) %>%
                  as.matrix(),
                lambda = ..6)
      ),
    model == "Regional Gtrends" ~
      pmap(
        list(data, arima_1, arima_2, arima_3, K, lambda),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$region_hits) %>%
                  as.matrix(),
                lambda = ..6)
      ),
    model == "FluVirus" ~
      pmap(
        list(data, arima_1, arima_2, arima_3, K, lambda),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$cum_h1per, ..1$cum_h3per) %>%
                  as.matrix(),
                lambda = ..6)
      ),
    model == "Backfill" ~
      pmap(
        list(data, arima_1, arima_2, arima_3, K, lambda),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$backfill) %>%
                  as.matrix(),
                lambda = ..6)
      ),
    model == "FluVirus & Backfill" ~
      pmap(
        list(data, arima_1, arima_2, arima_3, K, lambda),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$cum_h1per, ..1$cum_h3per,
                             ..1$backfill) %>%
                  as.matrix(),
                lambda = ..6)
      ),
    model == "National Gtrends & FluVirus" ~
      pmap(
        list(data, arima_1, arima_2, arima_3, K, lambda),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$hits, ..1$cum_h1per,
                             ..1$cum_h3per) %>%
                  as.matrix(),
                lambda = ..6)
      ),
    model == "Regional Gtrends & FluVirus" ~
      pmap(
        list(data, arima_1, arima_2, arima_3, K, lambda),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$region_hits, ..1$cum_h1per,
                             ..1$cum_h3per) %>%
                  as.matrix(),
                lambda = ..6)
      ),
    model == "National Gtrends & Backfill" ~
      pmap(
        list(data, arima_1, arima_2, arima_3, K, lambda),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$hits, ..1$backfill) %>%
                  as.matrix(),
                lambda = ..6)
      ),
    model == "Regional Gtrends & Backfill" ~
      pmap(
        list(data, arima_1, arima_2, arima_3, K, lambda),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$region_hits, ..1$backfill) %>%
                  as.matrix(),
                lambda = ..6)
      ),
    model == "National Gtrends, FluVirus, & Backfill" ~
      pmap(
        list(data, arima_1, arima_2, arima_3, K, lambda),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$hits, ..1$cum_h1per,
                             ..1$cum_h3per, ..1$backfill) %>%
                  as.matrix(),
                lambda = ..6)
      ),
    model == "Regional Gtrends, FluVirus, & Backfill" ~
      pmap(
        list(data, arima_1, arima_2, arima_3, K, lambda),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$region_hits, ..1$cum_h1per,
                             ..1$cum_h3per, ..1$backfill) %>%
                  as.matrix(),
                lambda = ..6)
      )
  )) %>%
  mutate(season = "2018/2019",
         prev_season = "2017/2018",
         week = EW,
         epiweek = epiweek)

# Create and save forecast files
springbok_path <- "State Forecasts/2018-2019/Dynamic Harmonic Model/"

set.seed(101085)
springbok_fit <- fits %>%
  mutate(
    pred_data = pmap(list(season, order_week, location, epiweek), 
                     ~ filter(state_flu_data_merge, year <= as.numeric(substr(..1, 6, 9)),
                              season != paste0(substr(..1, 6, 9), "/",
                                               as.numeric(substr(..1, 6, 9)) + 1),
                              season != ..1 | order_week %in% 40:..2,
                              location == ..3) %>%
                       left_join(select(state_ili_init_pub_list[[paste(..4)]],
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
    xreg = map(xreg,
               ~ as.matrix(.)),
    max_week = ifelse(season == "2014/2015", 53, 52),
    # Fit models
    pred_fit = pmap(list(pred_data, fit, xreg),
                    ~ Arima(..1$ILI, xreg = ..3, model = ..2)),
    # Create data frame of xreg terms for forecasting
    gtrend_forecast = pmap(
      list(pred_data, location, season, prev_season, order_week, max_week),
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
      list(pred_data, location, season, prev_season, order_week, max_week),
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
      list(data, order_week, max_week, season),
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
                        npaths = 1000,
                        challenge = "state_ili")
    )
  ) %>%
  select(location, week, max_week, pred_results)

springbok_pred <- springbok_fit %>%
  unnest() %>%
  select(-week, -max_week, -forecast_week) %>%
  bind_rows(generate_point_forecasts(.)) %>%
  select(location, target, type, unit, bin_start_incl, bin_end_notincl, value)

write_csv(springbok_pred,path = paste0(springbok_path, "/EW", EW_paste, ".csv"))

##### Cheetah #####

# List forecast files to be weighted
model_files <- list.files(path = "State Forecasts/2018-2019", recursive = TRUE, full.names = T)
model_files <- model_files[!grepl('ens', model_files)]

# List weighting schemes to be calculated
weight_files <- list.files(path = "State Weights")
weight_types <- sub(
  pattern = "-weights.csv", 
  replacement = "",
  sub(pattern = "Weights/",
      replacement = "",
      weight_files)
)

# Fetch ensemble with best performance for each state
load("Data/state_best_model_fits.Rdata")

final_sub <- tibble()

for(j in 1:length(weight_files)) {
  stacking_weights <- read.csv(paste0("State Weights/", weight_files[j]), 
                               stringsAsFactors=FALSE) %>%
    filter(season == "2018/2019")
  stacked_name <- sub(pattern = ".csv", replacement = "", weight_files[j])
  
  # If weights by week included, reset to MMWR week
  if (any(grepl("Model.Week", names(stacking_weights)))) {
    stacking_weights <- mutate(
      stacking_weights,
      `Model.Week` = week_reset(`Model.Week`, season),
      ew = paste0("EW", str_pad(`Model.Week`, 2, "left", 0))
    ) %>%
      select(-`Model.Week`) 
  }
  
  wt_subset <- select(stacking_weights, -season)
    
  dir.create(file.path("State Forecasts/2018-2019",
                       paste0("ens-", stacked_name)), 
             showWarnings = FALSE)
  
  ## identify the "EWXX-YYYY" combos for files given season
  # first_year <- 2018
  # first_year_season_weeks <- 40:52
  # week_names <- c(paste0("EW", first_year_season_weeks),
  #                 paste0("EW", str_pad(1:20, 2, "left", pad = 0)))
    
  this_week <- paste0("EW", EW_paste)
  
  if (any(grepl("ew", names(wt_subset)))) {
    wt_sub_subset <- filter(wt_subset, ew == this_week) %>%
      select(-ew)
  } else {
    wt_sub_subset <- wt_subset
  }
  
  ## stack models, save ensemble file
  files_to_stack <- model_files[grepl(this_week, model_files)]
  
  file_df <- data.frame(
    file = files_to_stack, 
    model_id = word(files_to_stack, start = 3, end = 3, sep = "/"),
    stringsAsFactors = FALSE)
  
  stacked_entry <- suppressMessages(
    stack_forecasts(file_df, wt_sub_subset, challenge = "state_ili")
  ) %>%
    select(location, target, type, unit, bin_start_incl, bin_end_notincl, value)
  
  # If this ensemble is the best for a state, add predictions for that state to final
  if(paste0("ens-", stacked_name) %in% scores_by_state$team) {
    final_sub <- bind_rows(
      final_sub,
      filter(stacked_entry, 
             location == scores_by_state$location[
               scores_by_state$team == paste0("ens-", stacked_name)])
    )
  }
  
  stacked_file_name <- paste0(
    "State Forecasts/2018-2019/ens-",
    stacked_name, "/", this_week, ".csv"
  )
  write.csv(stacked_entry, file=stacked_file_name, 
            row.names = FALSE, quote = FALSE)
  
}

# Write final submission CSV
write.csv(final_sub, 
          file = paste0("State Forecasts/2018-2019/ens-optimal-state/",
                        this_week, ".csv"),
          row.names = FALSE, quote = FALSE)


