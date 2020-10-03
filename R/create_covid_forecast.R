library(tidyverse)
library(FluSight)
library(forecast)

# Load data ------
ili_current <- readRDS("Data/ili_current.RDS")
ili_init_pub_list <- readRDS("Data/ili_init_pub_list.RDS")

covid_tracking <- readRDS('Data/covid_tracking.RDS') %>%
  filter(!location %in% c('District of Columbia', 'Puerto Rico')) %>%
  arrange(location, date) %>%
  group_by(location) %>%
  mutate(per_pos = coalesce(positive / (positive + negative), 0),
         new_pos = coalesce(positive - lag(positive), 0)) %>%
  ungroup() %>%
  left_join(select(ili_current, week, location, year, ILI),
            by = c('week', 'location', 'year'))

ggplot(filter(covid_tracking, !location %in% state.name), aes(x = date, y = per_pos)) +
  geom_line() +
  facet_wrap(~location)
ggplot(filter(covid_tracking, !location %in% state.name), aes(x = date, y = new_pos)) +
  geom_line() +
  facet_wrap(~location, scales = 'free')
ggplot(filter(covid_tracking, !location %in% state.name), aes(x = new_pos, y = ILI)) +
  geom_point() +
  facet_wrap(~location, scales = 'free')
ggplot(filter(covid_tracking, !location %in% state.name), aes(x = per_pos, y = ILI)) +
  geom_point() +
  facet_wrap(~location, scales = 'free')

ggplot(filter(ili_current, !location %in% state.name, year == 2019), aes(x = week, y = ILI)) +
  geom_line() +
  facet_wrap(~location)


covid_gtrend <- readRDS('Data/covid_gtrend.RDS') %>%
  bind_rows(
    readRDS('Data/covid_gtrend.RDS') %>%
      # Create HHS region proxies for Gtrend data
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
        TRUE ~ NA_character_
      )) %>%
      filter(!is.na(location))
  ) %>%
  filter(year > 2020 | (year == 2020 & week >= 4)) %>%
  # filter(location == 'US National') %>%
  # mutate(covid_hits = ts(covid_hits, frequency = 52, start = c(2020, 04)),
  #        coronavirus_hits = ts(coronavirus_hits, frequency = 52, start = c(2020, 04))) %>%
  left_join(select(ili_current, week, location, year, ILI),
            by = c('week', 'location', 'year'))

covid_fit <- auto.arima(covid_gtrend$covid_hits)
coronavirus_fit <- auto.arima(covid_gtrend$coronavirus_hits)
covid_test_fit <- auto.arima(covid_gtrend$covid_test_hits)

checkresiduals(covid_fit)
checkresiduals(coronavirus_fit)

autoplot(forecast(covid_fit))
autoplot(forecast(coronavirus_fit))
  # nest(data = -location) %>%
  # mutate(data = map(data, ~ mutate(., covid_hits = ts(covid_hits, frequency = 52, start = c(2020, 04)),
  #                                  coronavirus_hits = ts(coronavirus_hits, frequency = 52, start = c(2020, 04)))))
?ts
ggplot(filter(covid_gtrend, !location %in% state.name), aes(x = covid_hits, y = ILI)) +
  geom_point() +
  facet_wrap(~location)
ggplot(filter(covid_gtrend, !location %in% state.name), aes(x = coronavirus_hits, y = ILI)) +
  geom_point() +
  facet_wrap(~location)
ggplot(filter(covid_gtrend, !location %in% state.name), aes(x = covid_test_hits, y = ILI)) +
  geom_point() +
  facet_wrap(~location)


ggplot(filter(covid_gtrend, !location %in% state.name), aes(x = date, y = covid_hits)) +
  geom_point() +
  facet_wrap(~location)
?nest
?auto.arima()
