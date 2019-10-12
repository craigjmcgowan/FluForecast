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
  filter(!location %in% state.name)

ili_scored <- readRDS('Data/ili_init_pub_list.RDS')[c('201428', '201528', '201628', 
                                                      '201728', '201828', '201928')] %>%
  bind_rows() %>%
  filter(season != "2012/2013", !location %in% state.name) %>%
  bind_rows(filter(ili_current, 
                   season %in% c('2010/2011', '2011/2012', '2012/2013'))) 

ili_backfill <- readRDS('Data/ili_init_pub_list.RDS') %>%
  lapply(function(x) dplyr::filter(x, !location %in% state.name)) %>%
  lapply(add_scored_value) %>%
  bind_rows() %>%
  mutate(measure_week = week_inorder(week, season),
         week = measure_week + lag,
         print_year = as.numeric(substr(season, 1, 4))) %>%
  filter(lag <= 35, measure_week <= 73, week <= 73) %>%
  select(location, print_year, season, week, measure_week, backfill)


ili_backfill_densities <- ili_backfill %>%
  select(location, print_year, season, week, measure_week) %>%
  mutate(full_data = pmap(list(location, week, measure_week),
                          ~ filter(ili_backfill, location == ..1,
                                   week == ..2, measure_week == ..3) %>%
                            select(backfill)),
         sub_data = pmap(list(location, print_year, week, measure_week),
                         ~ filter(ili_backfill, location == ..1,
                                  print_year < max(..2, 2015), 
                                  week == ..3, measure_week == ..4) %>%
                           select(backfill)),
         data = ifelse(location != "US National" & week %in% c(40, 73),
                       full_data, sub_data),
         density_bw = map(data,
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
  left_join(filter(ili_backfill_densities, season == "2014/2015") %>%
               select(-season),
            by = c('location', 'week', 'measure_week'))

# Save densities
bind_rows(ili_backfill_densities, ili_fill_densities) %>%
  nest(backfill = c(measure_week, data, density_bw)) %>%
  saveRDS("Data/ili_backfill_densities.Rds")
