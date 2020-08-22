library(tidyverse)
library(lubridate)
library(MMWRweek)

source("R/utils.R")

### Covidtracking.com
covid_tracking_us <- read_csv('https://covidtracking.com/api/v1/us/daily.csv') %>%
  mutate(state = 'US')

covid_tracking_state <- read_csv('https://covidtracking.com/api/v1/states/daily.csv')

cols_keep <- names(covid_tracking_state)[names(covid_tracking_state) %in% names(covid_tracking_us)]

covid_tracking <- bind_rows(
  select(covid_tracking_us, all_of(cols_keep)),
  select(covid_tracking_state, all_of(cols_keep))
) %>%
  mutate(date = ymd(date),
         week = MMWRweek(date)$MMWRweek) %>%
  group_by(week, state) %>%
  arrange(date, .by_group = TRUE) %>%
  summarize(date = first(date),
            positive = coalesce(last(positive), 0),
            negative = coalesce(last(negative), 0),
            recovered = coalesce(last(recovered), 0),
            death = coalesce(last(death), 0),
            positiveIncrease = coalesce(mean(positiveIncrease, na.rm = T), 0),
            negativeIncrease = coalesce(mean(negativeIncrease, na.rm = T), 0)) %>%
  ungroup()

# Gtrends
covid_gtrend <- fetch_gtrend('US', 'covid')
coronavirus_gtrend <- fetch_gtrend('US', 'coronavirus')

ggplot(filter(coronavirus_gtrend, year == 2020), aes(x = date, y = hits)) +
  geom_line()
