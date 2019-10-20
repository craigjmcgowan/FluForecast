library(tidyverse)
library(zoo)
# library(cdcfluview)

source('R/utils.R')

percent_providers <- function(df, max_provider_df) {
  left_join(df, max_provider_df, by = c('location', 'season')) %>%
    mutate(per_report = num_providers / max_providers)
}

add_final_value <- function(df) {
  left_join(df, select(ili_current, season, location, week, final_ILI = ILI),
            by = c('season', 'location', 'week')) %>%
    mutate(backfill = final_ILI - ILI,
           per_backfill = final_ILI / ILI)
}

# Pull in data
ili_backfill <- readRDS('Data/ili_backfill.RDS') %>%
  filter(!location %in% state.name) %>%
  mutate(per_report = orig_num_providers / final_num_providers)

ili_current <- readRDS('Data/ili_current.RDS') %>%
  filter(!location %in% state.name) 

max_num_providers <- ili_current %>%
  group_by(season, location) %>%
  summarize(max_providers = max(num_providers))

ili_init_pub_list <- readRDS('Data/ili_init_pub_list.RDS') %>%
  lapply(function(x) dplyr::filter(x, !location %in% state.name)) %>%
  lapply(percent_providers, max_num_providers) %>%
  lapply(add_final_value)

# Look at relationships between backfill
ili_init_df <- bind_rows(ili_init_pub_list) %>%
  mutate(order_week = week_inorder(week, season))

ggplot(filter(ili_init_df, lag < 40), aes(x = week, y = per_backfill)) +
  scale_y_continuous(limits = c(0, 3)) +
  geom_point() +
  facet_wrap(~location)

ggplot(filter(ili_init_df, lag < 40), aes(x = week, y = backfill)) +
  geom_point() +
  facet_wrap(~location)

ggplot(filter(ili_init_df, lag < 40), aes(x = lag, y = per_backfill)) +
  geom_point() +
  facet_wrap(~location)

ggplot(ili_init_df, aes(x = per_report, y = backfill)) +
  geom_point() +
  facet_wrap(~location)

# Avg backfill by week/location/lag
avg_loc_wk_lag_backfill <- ili_init_df %>%
  filter(lag < 40, order_week >= 40 , order_week <= 74) %>%
  group_by(location, order_week, lag) %>%
  summarize(mean_backfill = mean(backfill),
            sd_backfill = sd(backfill),
            n = n())

ggplot(avg_loc_wk_lag_backfill, aes(x = n)) +
  geom_bar()

ggplot(avg_loc_wk_lag_backfill, aes(x = lag,  y = mean_backfill)) +
  geom_point() +
  facet_wrap(~location)

ggplot(avg_loc_wk_lag_backfill, aes(x = lag,  y = sd_backfill)) +
  geom_point() +
  facet_wrap(~location)

avg_loc_wk_lag_backfill %>%
  filter(n == 1) %>%
  head(20)

# Create backfill density similar to Delphi method
library(epiforecast)

full.dat = fetchEpidataFullDat("fluview", 'nat', "wili",
                               min.points.in.season = 52L,
                               first.week.of.season = 31L,
                               cache.file.prefix=sprintf("fluview_%s_fetch", 'nat'))
full.dat <- full.dat[names(full.dat)!="S2009"]
sim = twkde.markovian.sim(full.dat, baseline=2.1, max.n.sims=100)

new_sim = twkde.sim(full.dat, baseline=2.1, max.n.sims=100)

print(sim)
plot(sim)
