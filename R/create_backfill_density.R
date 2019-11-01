library(tidyverse)

source('R/utils.R')

# Helper functions
add_scored_value <- function(df) {
  inner_join(df, select(ili_scored, season, location, week, final_ILI = ILI),
            by = c('season', 'location', 'week')) %>%
    mutate(backfill = final_ILI - ILI,
           per_backfill = final_ILI / ILI)
}

fit_density <- function(x) {
  tryCatch(density(x, bw = "SJ")$bw,
           error=function(e) density(x)$bw)
}

# Load data and  create densities for each location, week, lag combination
ili_current <- readRDS('Data/ili_current.RDS') %>%
  mutate(release_date = as.Date(release_date)) %>%
  filter(season != "2019/2020")

ili_scored <- readRDS('Data/ili_init_pub_list.RDS')[c('201428', '201528', '201628', 
                                                      '201728', '201828', '201928')] %>%
  bind_rows() %>%
  filter(season != "2012/2013") %>%
  bind_rows(filter(ili_current, 
                   season %in% c('2010/2011', '2011/2012', '2012/2013'),
                   !location %in% state.name))

#####  Calculate backfill densities for national/regional #####
ili_backfill_nat_reg <- readRDS('Data/ili_init_pub_list.RDS') %>%
  lapply(function(x) dplyr::filter(x, !location %in% state.name)) %>%
  lapply(add_scored_value) %>%
  bind_rows() %>%
  mutate(measure_week = week_inorder(week, season),
         week = measure_week + lag,
         print_year = as.numeric(substr(season, 1, 4))) %>%
  filter(lag <= 35, measure_week <= 73, week <= 73) %>%
  select(location, print_year, season, week, measure_week, backfill)


ili_backfill_densities_nat_reg <- ili_backfill_nat_reg %>%
  # Add in rows for 2019/2020
  bind_rows(filter(ili_backfill_nat_reg, season == '2018/2019') %>%
              mutate(print_year = print_year + 1,
                     season = '2019/2020',
                     backfill = NA_real_)) %>%
  select(location, print_year, season, week, measure_week) %>%
  mutate(full_data = pmap(list(location, week, measure_week),
                          ~ filter(ili_backfill_nat_reg, location == ..1,
                                   week == ..2, measure_week == ..3) %>%
                            select(backfill)),
         sub_data = pmap(list(location, print_year, week, measure_week),
                         ~ filter(ili_backfill_nat_reg, location == ..1,
                                  print_year < max(..2, 2015), 
                                  week == ..3, measure_week == ..4) %>%
                           select(backfill)),
         data = ifelse(location != "US National" & week %in% c(40, 73),
                       full_data, sub_data),
         density_bw = map_dbl(data,
                          ~ fit_density(.$backfill))) %>%
  select(location, season, week, measure_week, data, density_bw) 

# Create entries for regional week 40s in early seasons
ili_fill_densities <- crossing(location = c("HHS Region 1", "HHS Region 2", "HHS Region 3",
                                        "HHS Region 4", "HHS Region 5", "HHS Region 6", 
                                        "HHS Region 7", "HHS Region 8",  "HHS Region 9",
                                        "HHS Region 10"),
                           season = c("2010/2011", "2011/2012", "2012/2013", "2013/2014"),
                           week = 40,
                           measure_week = 40) %>%
  left_join(filter(ili_backfill_densities_nat_reg, season == "2014/2015") %>%
               select(-season),
            by = c('location', 'week', 'measure_week'))

##### Calculate state backfill densities #####
ili_backfill_state <- readRDS('Data/ili_init_pub_list.RDS') %>%
  lapply(function(x) dplyr::filter(x, location %in% state.name)) %>%
  lapply(add_scored_value) %>%
  bind_rows() %>%
  mutate(measure_week = week_inorder(week, season),
         week = measure_week + lag,
         print_year = as.numeric(substr(season, 1, 4))) %>%
  filter(lag <= 35, measure_week <= 73, week <= 73) %>%
  select(location, print_year, season, week, measure_week, backfill)


# Identify locations with no backfill either year or one observation
ili_backfill_no_backfill <- ili_backfill_state %>%
  group_by(location, week, measure_week) %>%
  filter(near(min(backfill), max(backfill))) %>%
  select(location, week, measure_week) %>%
  distinct() %>%
  mutate(no_dens = 1)


ili_backfill_data_state <- ili_backfill_state %>%
  # Identify locations with one obs, no backfill, or identical values either year
  left_join(ili_backfill_no_backfill,
            by = c('location', 'week', 'measure_week')) %>%
  # Keep single location/week/lag combo - same estimate for all years
  select(location, week, measure_week, no_dens) %>%
  distinct() %>%
  mutate(data = pmap(list(location, week, measure_week),
                          ~ filter(ili_backfill_state, location == ..1,
                                   week == ..2, measure_week == ..3) %>%
                            select(backfill)))

ili_backfill_densities_state <- ili_backfill_data_state %>%
  filter(is.na(no_dens)) %>%
  mutate(density_bw = map_dbl(data, ~ fit_density(.$backfill))) %>%
  select(location, week, measure_week, density_bw) 

# Add BW equal to median for that location in for locations with no prior backfill or insufficient data
med_bw_by_state <- ili_backfill_densities_state %>%
  group_by(location) %>%
  summarize(med_bw = median(density_bw))


ili_backfill_densities_state_full <- ili_backfill_data_state %>%
  select(-no_dens) %>%
  left_join(ili_backfill_densities_state, 
            by = c('location', 'week', 'measure_week')) %>%
  left_join(med_bw_by_state, by = 'location') %>%
  mutate(density_bw = case_when(!is.na(density_bw) ~ density_bw,
                                !is.na(med_bw) ~ med_bw,
                                TRUE ~ median(ili_backfill_densities_state$density_bw)))  %>%
  select(-med_bw)

# Save densities
bind_rows(ili_backfill_densities_nat_reg, ili_fill_densities, 
          ili_backfill_densities_state_full) %>%
  nest(backfill = c(measure_week, data, density_bw)) %>%
  saveRDS("Data/ili_backfill_densities.Rds")



