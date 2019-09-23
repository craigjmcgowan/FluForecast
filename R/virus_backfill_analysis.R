library(tidyverse)
library(zoo)

source('R/utils.R')

# Explore backfill in virologic data

# Helper functions
subtype_read <- function(file) {
  read_csv(file) %>%
    mutate(epiweek = as.numeric(substr(basename(file), 3, 4)))
}

# Read in data
nat_files <- list.files('Data/WHO-NREVSS Data', "*Nat_PH*", full.names = TRUE)
reg_files <- list.files('Data/WHO-NREVSS Data', '*Reg_PH*', full.names = TRUE)

files <- c(nat_files, reg_files)

virus_data <- lapply(
  files,
  subtype_read
) %>%
  bind_rows()

# Create data frame of initial values
initial_values <- virus_data %>%
  filter(week >= 47 | week <= 18) %>%
  filter(week == epiweek) %>%
  select(region, year, week,
         total_specimens_init = total_specimens,
         a_2009_h1n1_init = a_2009_h1n1,
         a_h3_init = a_h3,
         a_subtyping_not_performed_init = a_subtyping_not_performed,
         b_init = b,
         bvic_init = bvic,
         byam_init = byam,
         h3n2v_init = h3n2v)

final_values <- virus_data %>%
  filter(epiweek == 37) %>%
  filter(week >= 47 | week <= 18) %>%
  mutate(order_week = week_inorder(week, year)) %>%
  group_by(region) %>%
  arrange(order_week, .by_group = TRUE) %>%
  mutate(h3_per = a_h3 / (a_2009_h1n1 + a_h3),
         h3_per = if_else(is.na(h3_per), 0.5, h3_per),
         h1_per = 1 - h3_per,
         h1_per_4wk = rollapply(h1_per, 4, mean, align = 'right', partial = TRUE),
         h3_per_4wk = 1 - h1_per_4wk,
         h1_per_6wk = rollapply(h1_per, 6, mean, align = 'right', partial = TRUE),
         h3_per_6wk = 1 - h1_per_6wk,
         h1_per_8wk = rollapply(h1_per, 8, mean, align = 'right', partial = TRUE),
         h3_per_8wk = 1 - h1_per_8wk,
         h1_per_cum = cummean(h1_per),
         h3_per_cum = 1 - h1_per_cum) %>%
  ungroup() %>%
  select(region, year, week,
         total_specimens_final = total_specimens,
         a_2009_h1n1_final = a_2009_h1n1,
         a_h3_final = a_h3,
         a_subtyping_not_performed_final = a_subtyping_not_performed,
         b_final = b,
         bvic_final = bvic,
         byam_final = byam,
         h3n2v_final = h3n2v,
         h1_per_final = h1_per,
         h3_per_final = h3_per,
         h1_per_4wk_final = h1_per_4wk,
         h3_per_4wk_final = h3_per_4wk,
         h1_per_6wk_final = h1_per_6wk,
         h3_per_6wk_final = h3_per_6wk,
         h1_per_8wk_final = h1_per_8wk,
         h3_per_8wk_final = h3_per_8wk,
         h1_per_cum_final = h1_per_cum,
         h3_per_cum_final = h3_per_cum)

# Join initial values back onto main data
compare_virus_data = virus_data %>%
  full_join(initial_values, by = c('region', 'year', 'week')) %>%
  full_join(final_values, by = c('region', 'year', 'week')) %>%
  filter(week >= 47 | week <= 18) %>%
  # Calculate absolute backfill
  mutate(specimens_diff = total_specimens_final - total_specimens,
         specimens_per = total_specimens / total_specimens_final,
         h3_per = a_h3 / (a_2009_h1n1 + a_h3),
         h3_per = if_else(is.na(h3_per), 0.5, h3_per),
         h1_per = 1 - h3_per,
         order_week = week_inorder(week, year),
         epiweek_in_order = week_inorder(epiweek, year)) %>%
  # Create rolling %
  group_by(region, epiweek_in_order) %>%
  arrange(order_week, .by_group = TRUE) %>%
  mutate(h1_per_4wk = rollapply(h1_per, 4, mean, align = 'right', partial = TRUE),
         h3_per_4wk = 1 - h1_per_4wk,
         h1_per_6wk = rollapply(h1_per, 6, mean, align = 'right', partial = TRUE),
         h3_per_6wk = 1 - h1_per_6wk,
         h1_per_8wk = rollapply(h1_per, 8, mean, align = 'right', partial = TRUE),
         h3_per_8wk = 1 - h1_per_8wk,
         h1_per_cum = cummean(h1_per),
         h3_per_cum = 1 - h1_per_cum) %>%
  ungroup() %>%
  # Calculate differences in rolling averages
  mutate(h1_4wk_diff = h1_per_4wk - h1_per_4wk_final,
         h3_4wk_diff = h3_per_4wk - h3_per_4wk_final,
         h1_6wk_diff = h1_per_6wk - h1_per_6wk_final,
         h3_6wk_diff = h3_per_6wk - h3_per_6wk_final,
         h1_8wk_diff = h1_per_8wk - h1_per_8wk_final,
         h3_8wk_diff = h3_per_8wk - h3_per_8wk_final,
         h1_cum_diff = h1_per_cum - h1_per_cum_final,
         h3_cum_diff = h3_per_cum - h3_per_cum_final)
  
  
# Exploration
compare_virus_data %>%
  filter(week == epiweek) %>%
  ggplot(aes(x = wk_date, y = specimens_per)) +
  geom_line() +
  facet_wrap(~ region, ncol = 3)

filter(compare_virus_data, week == epiweek) %>%
  select(wk_date, region, h1_per, h1_per_final) %>%
  gather(key = 'Measurement', value = 'value', h1_per, h1_per_final) %>%
  mutate(Measurement = if_else(Measurement == "h1_per", "Initial", "Final")) %>%
  ggplot(aes(x = wk_date, y = value, color = Measurement)) +
  geom_line() +
  facet_wrap(~ region, ncol = 3)

# Rolling backfill
filter(compare_virus_data, week == epiweek) %>%
  select(wk_date, region, h1_4wk_diff, h1_6wk_diff, h1_8wk_diff, h1_cum_diff) %>%
  gather(key = "window", value = 'value', -wk_date, -region) %>%
  ggplot(aes(x = wk_date, y = value, color = window)) +
  geom_hline(yintercept = 0) +
  geom_line() +
  facet_wrap(~ region, ncol = 3)

filter(compare_virus_data, week == epiweek) %>%
  select(region, wk_date, h1_per_4wk, h1_per_4wk_final) %>%
  gather(key = 'measurement', value = 'value', h1_per_4wk, h1_per_4wk_final) %>%
  ggplot(aes(x = wk_date, y = value, color = measurement)) +
  geom_line() +
  facet_wrap(~ region)

filter(compare_virus_data, week == epiweek) %>%
  select(region, wk_date, h1_per_6wk, h1_per_6wk_final) %>%
  gather(key = 'measurement', value = 'value', h1_per_6wk, h1_per_6wk_final) %>%
  ggplot(aes(x = wk_date, y = value, color = measurement)) +
  geom_line() +
  facet_wrap(~ region)

filter(compare_virus_data, week == epiweek) %>%
  select(region, wk_date, h1_per_8wk, h1_per_8wk_final) %>%
  gather(key = 'measurement', value = 'value', h1_per_8wk, h1_per_8wk_final) %>%
  ggplot(aes(x = wk_date, y = value, color = measurement)) +
  geom_line() +
  facet_wrap(~ region)

filter(compare_virus_data, week == epiweek) %>%
  select(region, wk_date, h1_per_cum, h1_per_cum_final) %>%
  gather(key = 'measurement', value = 'value', h1_per_cum, h1_per_cum_final) %>%
  ggplot(aes(x = wk_date, y = value, color = measurement)) +
  geom_line() +
  facet_wrap(~ region)

# MAE
mae <- compare_virus_data %>%
  filter(week == epiweek) %>%
  group_by(region) %>%
  summarize(`4` = mean(abs(h1_4wk_diff)),
            `6` = mean(abs(h1_6wk_diff)),
            `8` = mean(abs(h1_8wk_diff)),
            `10` = mean(abs(h1_cum_diff)))

gather(mae, key = 'delay', value = 'value', -region) %>%
  mutate(delay = as.numeric(delay)) %>%
  ggplot(aes(x = delay, y = value, color = region)) +
  geom_line()


region_4 <- filter(compare_virus_data, region == "Region 4")
