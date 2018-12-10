library(tidyverse)
library(forecast)

load("Data/CV_ARIMA_terms.Rdata")
load("Data/CV_Fourier_terms.Rdata")
load("Data/CV_Transform_terms.Rdata")
load("Data/CV_covar_terms.Rdata")

load("Data/ili.Rdata")
load("Data/virologic.Rdata")
load("Data/Gtrends.Rdata")
load("Data/state_gtrend.Rdata")

flu_data_merge <- select(ili_current, epiweek, ILI, year, week, season, location) %>%
  inner_join(select(virologic_combined, location, season, year, week, cum_h1per,
                   cum_h3per, cum_bper),
            by = c("location", "season", "year", "week")) %>%
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
  right_join(gtrend_US_flu_merge,
            by = c("season", "year", "week")) %>%
  full_join(select(ili_backfill, -orig_ILI, -final_ILI),
            by = c("location", "season", "year", "week")) %>%
  # Remove 2008/2009 and 2009/2010 seasons due to pandemic activity
  filter(!season %in% c("2008/2009", "2009/2010")) %>%
  # Remove week 33 in 2014 so all seasons have 52 weeks - minimal activity
  filter(!(year == 2014 & week == 33)) %>%
  mutate(order_week = case_when(
    week < 40 & season == "2014/2015" ~ week + 53,
    week < 40 ~ week + 52,
    TRUE ~ week
  )) %>%
  # Replace missing backfill data with random draw from same season/week observed values
  mutate(
    sim_backfill = map2_dbl(location, week,
                        ~ sample(ili_backfill$backfill[ili_backfill$location == .x &
                                                         ili_backfill$week == .y], 1)),
    backfill = case_when(
      !is.na(backfill) ~ backfill,
      TRUE ~ sim_backfill
      )
    ) %>%
  select(-sim_backfill)



temp <- flu_data_merge %>%
  mutate(log_ILI = ifelse(ILI == 0, log(0.1), log(ILI)),
         week = ifelse(year == 2014 & week > 33, week - 1, week),
         year_frac = year + week/52,
         cos1 = cos(2*year_frac*pi/1),
         sin1 = sin(2*year_frac*pi/1),
         cos2 = cos(2*year_frac*pi/0.5),
         sin2 = sin(2*year_frac*pi/0.5)) %>%
  nest(-location) %>%
  filter(location == "US National") %>%
  mutate(data = map(data,
                    ~ mutate(., log_ILI_ts = ts(log_ILI, frequency = 52,
                                                start = c(2004, 40)),
                             fourier_1 = fourier(log_ILI_ts, K = 1)[, 1],
                             fourier_2 = fourier(log_ILI_ts, K = 1)[, 2])),
         xreg = map(data,
                    ~ select(., fourier_1, fourier_2, hits, cum_h1per,
                             cum_h3per)),
         fourier_mod_1 = map(data,
                             ~ glm(log_ILI ~ cos1 + sin1, data = .)),
         fourier_mod_2 = map(data,
                             ~ glm(log_ILI ~ cos1 + sin1 + cos2 + sin2,
                                  data = .)),
         # nat_sarima_mod = map(data,
         #                   ~ Arima(.x$log_ILI_ts, order = c(2, 0, 2),
         #                           seasonal = c(2, 1, 0))),
         nat_arima_mod = map2(data, xreg,
                             ~ Arima(.x$log_ILI_ts, order = c(0, 1, 3),
                                     xreg = .y)),
         log_pred_four_1 = pmap(list(fourier_mod_1, data),
                            ~ predict(..1, ..2, type='response')),
         log_pred_four_2 = pmap(list(fourier_mod_2, data),
                            ~ predict(..1, ..2, type='response')),
         # log_pred_sarima = pmap(list(nat_sarima_mod),
         #                        ~ fitted(..1)),
         log_pred_arima = pmap(list(nat_arima_mod, xreg),
                               ~ fitted(..1, xreg = ..2))) %>%
  select(-fourier_mod_1, -fourier_mod_2, -nat_arima_mod) %>% #, -nat_sarima_mod) %>%
  unnest() %>%
  mutate(pred_four_1 = exp(log_pred_four_1),
         pred_four_2 = exp(log_pred_four_2),
         # pred_nat_s = exp(as.numeric(log_pred_sarima)),
         pred_nat_arima = exp(as.numeric(log_pred_arima)))

ggplot(temp) +
  geom_point(aes(x = date, y = ILI)) +
  geom_line(aes(x = date, y = pred_four_1), color = "red", 
            size = 2, alpha = 0.5) +
  geom_line(aes(x = date, y = pred_four_2), color = "blue",
            size = 2, alpha = 0.5) +
  facet_wrap(~ location) +
  theme_minimal()

p <- ggplot(filter(temp, location == "US National")) +
  geom_line(aes(x = date, y = ILI)) +
  facet_wrap(~ season, scales = "free_x") + 
  theme_minimal()

p

p + geom_line(aes(x = date, y = pred_four_1), color = "red")

p + geom_line(aes(x = date, y = pred_four_2), color = "red")

p + geom_line(aes(x = date, y = pred_four_2), color = "red") +
  geom_line(aes(x = date, y = pred_nat_arima), color = "blue")
