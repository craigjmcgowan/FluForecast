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
EW <- 01
epiweek <- 201901

##### Update data #####
source("R/read_data.R")
source("R/gtrends_pull_data.R")
source("R/save_who_nrevss.R")

set.seed(4321) # For reproducibility of backfill sample

flu_data_merge <- select(ili_current, epiweek, ILI, year, week, season, location) %>%
  inner_join(select(virologic_combined, location, season, year, week, cum_h1per,
                   cum_h3per),
            by = c("location", "season", "year", "week")) %>%
  # Join backfill and replace missing backfill data with random draw from 
  # same season/week observed values
  full_join(select(ili_backfill, -orig_ILI, -final_ILI),
            by = c("location", "season", "year", "week")) %>%
  mutate(
    sim_backfill = map2_dbl(location, week,
                            ~ sample(ili_backfill$backfill[ili_backfill$location == .x &
                                                             ili_backfill$week == .y], 1)),
    backfill = case_when(
      !is.na(backfill) ~ backfill,
      TRUE ~ sim_backfill
    )
  ) %>%
  select(-sim_backfill) %>%
  # Add Google Trends data by location
  right_join(bind_rows(mutate(gtrend_US_flu_merge, location = "US National"),
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
  right_join(gtrend_US_flu_merge,
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
  

##### Kudu #####
vir_ssn_per <- virologic_combined %>%
  group_by(season, location) %>%
  summarize(h1sum = sum(a_2009_h1n1, na.rm = T) + sum(a_h1, na.rm = T),
            h3sum = sum(a_h3, na.rm = T),
            h1per = h1sum / (h1sum + h3sum),
            h3per = h3sum / (h1sum + h3sum)) %>%
  select(-h1sum, -h3sum)

onsets <- ili_current %>%
  filter(year >= 2007, !season %in% c("2006/2007", "2009/2010"),
         (week >= 40 | week <= 20)) %>%
  group_by(location, season) %>%
  do(create_onset(., region = .$location[1], 
                  year = as.numeric(substr(.$season[1], 1, 4)))) %>%
  ungroup() %>%
  left_join(vir_ssn_per, by = c("location", "season"))

# Weight probability of no onset by virus % in seasons w/o onset
prob_no_onset <- onsets %>%
  group_by(location) %>%
  arrange(season, .by_group = T) %>%
  mutate(prob_no_onset = lag(cumsum(bin_start_incl == "none") / row_number(), default = 0)) %>%
  group_by(location, isna = (bin_start_incl != "none")) %>%
  mutate(h1_prob_no_onset = ifelse(isna, NA, cummean(h1per)),
         h3_prob_no_onset = ifelse(isna, NA, cummean(h3per))) %>%
  group_by(location) %>%
  mutate(h1_prob_no_onset = na.locf(h1_prob_no_onset, na.rm = FALSE),
         h3_prob_no_onset = na.locf(h3_prob_no_onset, na.rm = FALSE)) %>%
  mutate(h1_prob_no_onset = ifelse(is.na(h1_prob_no_onset), 0, 
                                   h1_prob_no_onset * prob_no_onset),
         h3_prob_no_onset = ifelse(is.na(h3_prob_no_onset), 0, 
                                   h3_prob_no_onset * prob_no_onset)) %>%
  ungroup() %>%
  select(location, season, prob_no_onset, h1_prob_no_onset, h3_prob_no_onset)

# Create forecasts 
kudu_ili <- ili_current %>%
  filter(year <= 2018, season != "2018/2019")

kudu_virologic <- virologic_combined %>%
  filter(season == "2018/2019") %>%
  # Create variable to keep weeks in order
  mutate(order_week = week_inorder(week, season)) %>%
  select(year,location, week, h1per = h1per_of_a, 
         h3per = h3per_of_a, order_week) %>%
  select(location, order_week, h1per, h3per)

# Create directory to store forecasts
dir.create("Forecasts/2018-2019/Subtype Historical Average/",
           showWarnings = FALSE)

kudu_flusight_path <- 
  "../cdc-flusight-ensemble/model-forecasts/real-time-component-models/Protea_Kudu" 

# Create target densities and functions
subtype_densities_1819 <- create_subtype_densities(
  kudu_ili, vir_ssn_per
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
    prob_no_onset = filter(prob_no_onset, season == "2018/2019")
  ) %>%
  select(location, target, type, unit, bin_start_incl, bin_end_notincl, value)
  
EW_paste <- str_pad(EW, 2, pad = "0")
  
write_csv(kudu_pred,
          path = paste0("Forecasts/2018-2019/Subtype Historical Average/EW", EW_paste, ".csv"))
  

write_csv(kudu_pred, 
          path = paste0(kudu_flusight_path, "/EW", EW_paste,
                        "-", substr(epiweek, 1, 4), "-Protea_Kudu.csv"))


##### Springbok #####
load("Data/CV_Transform_terms.Rdata")
load("data/CV_Fourier_terms.Rdata")
load("Data/CV_ARIMA_terms.Rdata")
load("Data/CV_covar_terms.Rdata")

fits <- filter(flu_data_merge, year <= 2018,
               season != "2018/2019") %>%
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
  mutate(season = "2018/2019",
         prev_season = "2017/2018",
         week = EW,
         order_week = if_else(EW < 40, EW + 52, EW),
         epiweek = epiweek)

# Create and save forecast files
springbok_path <- "Forecasts/2018-2019/Dynamic Harmonic Model/"

springbok_flusight_path <- 
  "../cdc-flusight-ensemble/model-forecasts/real-time-component-models/Protea_Springbok" 

cdc_path <- "CDC Submissions/2018-2019"

set.seed(4321)
springbok_fit <- fits %>%
  mutate(
    pred_data = pmap(list(season, order_week, location, epiweek), 
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
      list(pred_data, location, season, prev_season, order_week, max_week),
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
      list(pred_data, location, season, prev_season, order_week, max_week),
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
      list(data, order_week, max_week, season),
      function(data, order_week, max_week, season) {
        out <- numeric()
        for(i in order_week:(max_week + 24)) {
          
          temp <- data[data$order_week == i &
                           !data$season %in% c("2004/2005", "2005/2006", "2006/2007",
                                               "2008/2009"), ]
          
          out <- c(out, predict(lm(backfill ~ year, 
                                   data = temp),
                                data.frame(year = case_when(
                                  i <= max_week ~ as.numeric(substr(season, 1, 4)),
                                  TRUE ~ as.numeric(substr(season, 6, 9))
                                  ))))
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
  select(location, week, max_week, pred_results)

springbok_pred <- springbok_fit %>%
  unnest() %>%
  select(-week, -max_week, -forecast_week) %>%
  bind_rows(generate_point_forecasts(.)) %>%
  select(location, target, type, unit, bin_start_incl, bin_end_notincl, value)

write_csv(springbok_pred,path = paste0(springbok_path, "/EW", EW_paste, ".csv"))

write_csv(springbok_pred, 
          path = paste0(springbok_flusight_path, "/EW", EW_paste,
                        "-", substr(epiweek, 1, 4), "-Protea_Springbok.csv"))

write_csv(springbok_pred, 
          path = paste0(cdc_path, "/EW", EW_paste,
                        "-Protea_Springbok-", Sys.Date(), ".csv"))

##### Cheetah #####

cheetah_flusight_path <- 
  "../cdc-flusight-ensemble/model-forecasts/real-time-component-models/Protea_Cheetah"

# List forecast files to be weighted
model_files <- list.files(path = "Forecasts/2018-2019", recursive = TRUE, full.names = T)
model_files <- model_files[!grepl('ens', model_files)]

# List weighting schemes to be calculated
weight_files <- list.files(path = "Weights")
weight_types <- sub(
  pattern = "-weights.csv", 
  replacement = "",
  sub(pattern = "Weights/",
      replacement = "",
      weight_files)
)


for(j in 1:length(weight_files)) {
  stacking_weights <- read.csv(paste0("weights/", weight_files[j]), 
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
    
  dir.create(file.path("Forecasts/2018-2019",
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
  
  stacked_entry <- stack_forecasts(file_df, wt_sub_subset) %>%
    select(location, target, type, unit, bin_start_incl, bin_end_notincl, value)
  stacked_file_name <- paste0(
    "Forecasts/2018-2019/ens-",
    stacked_name, "/", this_week, ".csv"
  )
  write.csv(stacked_entry, file=stacked_file_name, 
            row.names = FALSE, quote = FALSE)
  
  # Save Cheetah model to FSN path as well
  if(stacked_name == "month-target-type-based-weights") {
    write.csv(stacked_entry, 
              file = paste0(cheetah_flusight_path, "/", this_week,
                            "-", substr(epiweek, 1, 4), "-Protea_Cheetah.csv"),
              row.names = FALSE, quote = FALSE)
    
    write.csv(stacked_entry, 
              file = paste0(cdc_path, "/", this_week,
                            "-Protea_Cheetah-", Sys.Date(), ".csv"),
              row.names = FALSE, quote = FALSE)
  }
}

# Recreate README with updated forecasts
rmarkdown::render("README.Rmd")
