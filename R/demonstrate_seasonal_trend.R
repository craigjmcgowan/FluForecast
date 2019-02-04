library(tidyverse)
library(MMWRweek)
# library(forecast)

# load("Data/CV_ARIMA_terms.Rdata")
# load("Data/CV_Fourier_terms.Rdata")
# load("Data/CV_Transform_terms.Rdata")
# load("Data/CV_covar_terms.Rdata")

load("Data/ili.Rdata")
# load("Data/virologic.Rdata")
# load("Data/Gtrends.Rdata")
# load("Data/state_gtrend.Rdata")

# Read training data CSVs and extract 1 wk ahead point forecasts
training_points <- tibble()

for(this_season in c("2010-2011", "2011-2012", "2012-2013", "2013-2014",
                     "2014-2015", "2015-2016", "2016-2017", "2017-2018")) {
  
  training_files <- list.files(path = file.path("Forecasts", this_season,
                                                "ens-month-target-type-region-based-weights"),
                               full.names = TRUE)
  
  
  for(j in seq_along(training_files)) {
    training_points <- read_csv(training_files[j]) %>%
      filter(type == "Point", target == "1 wk ahead") %>%
      mutate(week = as.numeric(substr(training_files[j], 66, 67)) + 1,
             season = this_season,
             order_week = case_when(
               this_season == "2014-2015" & week < 40 ~ week + 53,
               week < 40 ~ week + 52,
               TRUE ~ week
             ),
             year = case_when(
               week < 40 ~ as.numeric(substr(season, 6, 9)),
               TRUE ~ as.numeric(substr(season, 1, 4))
             )) %>%
      select(season, location, year, week, order_week, value) %>%
      bind_rows(training_points, .)
  }
}


training_plot_data <- ili_current %>%
  select(season, location, week, ILI) %>%
  mutate(season = str_replace(season, "/", "-"),
         order_week = case_when(
           season == "2014-2015" & week < 40 ~ week + 53,
           week < 40 ~ week + 52,
           TRUE ~ as.numeric(week)
         )) %>%
  filter(season %in% c("2010-2011", "2011-2012", "2012-2013", "2013-2014",
                       "2014-2015", "2015-2016", "2016-2017", "2017-2018")) %>%
  inner_join(training_points, by = c("season", "location", "week", "order_week")) %>%
  mutate(date = MMWRweek2Date(year, week)) %>%
  gather(key = "data_source", value = "value", ILI, value) %>%
  mutate(data_source = case_when(
    data_source == "ILI" ~ "Observed activity",
    data_source == "value" ~ "1 wk predicted activity"
  ))



ggplot(data = filter(training_plot_data, location == "US National")) +
  facet_wrap(~ season, scales = "free_x") +
  geom_line(aes(x = date, y = value, color = data_source)) +
  scale_color_manual(values = c("red", "black")) +
  labs(x = "Date", y = "Percent outpatient visits due to ILI") + 
  theme_minimal()


# Several 4 wk plots with prediction intervals and observed data
ili_1819 <- ili_current %>%
  filter(season == "2017/2018", location == "US National") %>%
  select(season, location, week, ILI)

wk42 <- read_csv("CDC Submissions/2018-2019/EW42-Protea_Cheetah-2018-10-29.csv") %>%
  filter(location == "US National", target %in% c("1 wk ahead", "2 wk ahead",
                                                  "3 wk ahead", "4 wk ahead")) %>%
  mutate(week = 42,
         season = this_season,
         order_week = case_when(
           this_season == "2014-2015" & week < 40 ~ week + 53,
           week < 40 ~ week + 52,
           TRUE ~ week
         ),
         year = case_when(
           week < 40 ~ as.numeric(substr(season, 6, 9)),
           TRUE ~ as.numeric(substr(season, 1, 4))
         )) %>%
  select(season, location, year, week, order_week, value)


# Dump regional plot here for now
# ggplot(filter(plot_points, forecast_date == as.Date("2018-12-09"),
#               location != "US National")) +
#   # Add shading for 80% CI for average
#   geom_ribbon(data = filter(bounds, forecast_date == as.Date("2018-12-09"),
#                             location != "US National"),
#               aes(x = date, ymin = lower, ymax = upper, 
#                   fill = "80% prediction interval"),
#               alpha = 0.3) +
#   # Add point prediction lines for each team
#   geom_line(aes(date, value, color = "Predicted values")) +
#   # Add observed values
#   geom_line(data = filter(plot_model_truth, forecast_date == as.Date("2018-12-09"),
#                           location != "US National"),
#             aes(date, value, linetype = future_ILI), color = "black") +
#   # Make pretty
#   theme_minimal() +
#   theme(legend.position = "bottom") +
#   labs(x = "Date", y = "Percent outpatient visits due to ILI",
#        linetype = "", color = "", fill = "") +
#   guides(linetype = guide_legend(order = 1),
#          color = guide_legend(order = 2),
#          fill = guide_legend(order = 3)) +
#   facet_wrap( ~ location, scales = "free_y", nrow = 2)
