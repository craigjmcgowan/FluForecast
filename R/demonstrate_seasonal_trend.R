library(tidyverse)
library(forecast)


load("Data/CV_ARIMA_terms.Rdata")
load("Data/CV_Fourier_terms.Rdata")
load("Data/CV_Transform_terms.Rdata")
load("Data/CV_covar_terms.Rdata")

load("Data/truth_and_data.Rdata")

temp <- flu_data_merge %>%
  mutate(log_ILI = ifelse(ILI == 0, log(0.1), log(ILI)),
         year_frac = year + week/52.17,
         cos1 = cos(2*year_frac*pi/1),
         sin1 = sin(2*year_frac*pi/1),
         cos2 = cos(2*year_frac*pi/0.5),
         sin2 = sin(2*year_frac*pi/0.5)) %>%
  nest(-location) %>%
  mutate(fourier_mod_1 = map(data,
                             ~ glm(log_ILI ~ cos1 + sin1, data = .)),
         fourier_mod_2 = map(data,
                             ~ glm(log_ILI ~ cos1 + sin1 + cos2 + sin2,
                                  data = .)),
         nat_mod = map(data,
                           ~ glm(log_ILI ~ cos1 + sin1 + hits + cum_h1per +
                                   cum_h3per, data = .)),
         nat_mod_back = map(data,
                       ~ glm(log_ILI ~ cos1 + sin1 + hits + cum_h1per +
                               cum_h3per + backfill, data = .)),
         log_pred_four_1 = pmap(list(fourier_mod_1, data),
                            ~ predict(..1, ..2, type='response')),
         log_pred_four_2 = pmap(list(fourier_mod_2, data),
                            ~ predict(..1, ..2, type='response')),
         log_pred_nat = pmap(list(nat_mod, data),
                                ~ predict(..1, ..2, type='response')),
         log_pred_nat_back = pmap(list(nat_mod_back, data),
                             ~ predict(..1, ..2, type='response'))) %>%
  select(-fourier_mod_1, -fourier_mod_2, -nat_mod, -nat_mod_back) %>%
  unnest() %>%
  mutate(pred_four_1 = exp(log_pred_four_1),
         pred_four_2 = exp(log_pred_four_2),
         pred_nat = exp(log_pred_nat),
         pred_nat_back = exp(log_pred_nat_back))

ggplot(temp) +
  geom_point(aes(x = date, y = ILI)) +
  geom_line(aes(x = date, y = pred_four_1), color = "red", 
            size = 2, alpha = 0.5) +
  geom_line(aes(x = date, y = pred_four_2), color = "blue",
            size = 2, alpha = 0.5) +
  facet_wrap(~ location) +
  theme_minimal()

ggplot(filter(temp, location == "US National")) +
  geom_line(aes(x = date, y = ILI)) +
  # geom_line(aes(x = date, y = pred_nat), color = "red", 
  #           size = 2, alpha = 0.5) +
  geom_line(aes(x = date, y = pred_nat_back), color = "blue", 
            size = 2, alpha = 0.5) +
  theme_minimal()

?geom_smooth

?predict.glm
?map2
