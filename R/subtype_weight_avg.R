# Generate subtype-weighted historical average forecast
library(tidyverse)
library(FluSight)
library(zoo)

# Load data -------
ili_current <- readRDS("Data/ili_current.RDS")
virologic_combined <- readRDS("Data/virologic.RDS")
pseudo_onsets <- readRDS("Data/pseudo_onsets.RDS")

# Load functions 
source("R/create_subtype_forecast.R")
source("R/utils.R")

kudu_path <- "../cdc-flusight-ensemble/model-forecasts/component-models/Protea_Kudu" 

# Calculate cumulative virologic status across seasons -----
vir_ssn_per <- virologic_combined %>%
  filter(!location %in% state.name) %>%
  mutate(order_week = week_inorder(week, season)) %>%
  group_by(season, location) %>%
  arrange(order_week, .by_group = TRUE) %>% 
  summarize(h1per = last(cum_h1per),
            h3per = last(cum_h3per),
            bper = last(cum_bper)) %>%
  ungroup()

vir_wk_per <- virologic_combined %>%
  filter(!location %in% state.name) %>%
  mutate(wk_date = MMWRweek::MMWRweek2Date(year, week),
         order_week = week_inorder(week, season)) %>%
  select(location, season, order_week, 
         h1per = cum_h1per_6wk, 
         h3per = cum_h3per_6wk,
         bper = cum_bper_6wk)

# Create data.frame of probability of no onset based on prior years -----
onsets <- ili_current %>%
  filter(year >= 2007, !season %in% c("2006/2007", "2009/2010"),
         (week >= 40 | week <= 20), !location %in% state.name) %>%
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
         h3_prob_no_onset = ifelse(isna, NA, cummean(h3per)),
         b_prob_no_onset = ifelse(isna, NA, cummean(bper))) %>%
  group_by(location) %>%
  mutate(h1_prob_no_onset = na.locf(h1_prob_no_onset, na.rm = FALSE),
         h3_prob_no_onset = na.locf(h3_prob_no_onset, na.rm = FALSE),
         b_prob_no_onset = na.locf(b_prob_no_onset, na.rm = FALSE)) %>%
  mutate(h1_prob_no_onset = ifelse(is.na(h1_prob_no_onset), 0, 
                                   h1_prob_no_onset * prob_no_onset),
         h3_prob_no_onset = ifelse(is.na(h3_prob_no_onset), 0, 
                                   h3_prob_no_onset * prob_no_onset),
         b_prob_no_onset = ifelse(is.na(b_prob_no_onset), 0, 
                                   b_prob_no_onset * prob_no_onset)) %>%
  ungroup() %>%
  select(location, season, prob_no_onset, h1_prob_no_onset, h3_prob_no_onset,
         b_prob_no_onset)

# Create forecasts for 2010/2011 ------
train_ili_1011 <- ili_current %>%
  filter(year <= 2010, season != "2010/2011",
         !location %in% state.name)

virologic_1011 <- virologic_combined %>%
  filter(season == "2010/2011", !location %in% state.name) %>%
  # Create variable to keep weeks in order
  mutate(order_week = week_inorder(week, season)) %>%
  select(year, location, week, order_week, 
         h1per_wk = cum_h1per_6wk, h3per_wk = cum_h3per_6wk,
         bper_wk = cum_bper_6wk, h1per_ssn = cum_h1per,
         h3per_ssn = cum_h3per, bper_ssn = cum_bper)

# Create directory to store forecasts
dir.create("Forecasts/Training/2010-2011/Subtype Historical Average/",
           showWarnings = FALSE)

# Create target densities and functions
subtype_densities_1011 <- create_subtype_densities(
  ili_df = train_ili_1011, vir_ssn_per = vir_ssn_per, vir_wk_per = vir_wk_per, 
  pseudo_onsets = pseudo_onsets
)
subtype_functions_1011 <- modify_depth(
  subtype_densities_1011, 3,
  function(dens) approxfun(dens$x, dens$y, rule = 2)
)

for(i in 40:72) {
  temp <- create_subtype_forecast(
    functions = subtype_functions_1011,
    virologic = virologic_1011,
    pub_week = i,
    season = "2010/2011",
    prob_no_onset = filter(prob_no_onset, season == "2010/2011")
  ) 
  
  temp$Value <- format(temp$Value, scientific = FALSE)
  
  j <- str_pad(ifelse(i > 52, i - 52, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("Forecasts/Training/2010-2011/Subtype Historical Average/EW", j, ".csv"))
  
  year <- ifelse(i > 52, 2011, 2010)
  write_csv(temp,
            path = paste0(kudu_path, "/EW", j, "-", year, "-Protea_Kudu.csv"))
  
}

# Create forecasts for 2011/2012 ------
train_ili_1112 <- ili_current %>%
  filter(year <= 2011, season != "2011/2012",
         !location %in% state.name)

virologic_1112 <- virologic_combined %>%
  filter(season == "2011/2012", !location %in% state.name) %>%
  # Create variable to keep weeks in order
  mutate(order_week = week_inorder(week, season)) %>%
  select(year, location, week, order_week, 
         h1per_wk = cum_h1per_6wk, h3per_wk = cum_h3per_6wk,
         bper_wk = cum_bper_6wk, h1per_ssn = cum_h1per,
         h3per_ssn = cum_h3per, bper_ssn = cum_bper)

# Create directory to store forecasts
dir.create("Forecasts/Training/2011-2012/Subtype Historical Average/",
           showWarnings = FALSE)

# Create target densities and functions
subtype_densities_1112 <- create_subtype_densities(
  train_ili_1112, vir_ssn_per = vir_ssn_per, vir_wk_per = vir_wk_per, 
  pseudo_onsets = pseudo_onsets
)

subtype_functions_1112 <- modify_depth(
  subtype_densities_1112, 3,
  function(dens) approxfun(dens$x, dens$y, rule = 2)
)

for(i in 40:72) {
  temp <- create_subtype_forecast(
    functions = subtype_functions_1112,
    virologic = virologic_1112,
    pub_week = i,
    season = "2011/2012",
    prob_no_onset = filter(prob_no_onset, season == "2011/2012")
  )
  
  temp$Value <- format(temp$Value, scientific = FALSE)
  
  j <- str_pad(ifelse(i > 52, i - 52, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("Forecasts/Training/2011-2012/Subtype Historical Average/EW", j, ".csv"))
  
  year <- ifelse(i > 52, 2012, 2011)
  write_csv(temp,
            path = paste0(kudu_path, "/EW", j, "-", year, "-Protea_Kudu.csv"))
  
}

# Create forecasts for 2012/2013 ------
train_ili_1213 <- ili_current %>%
  filter(year <= 2012, season != "2012/2013",
         !location %in% state.name)

virologic_1213 <- virologic_combined %>%
  filter(season == "2012/2013", !location %in% state.name) %>%
  # Create variable to keep weeks in order
  mutate(order_week = week_inorder(week, season)) %>%
  select(year, location, week, order_week, 
         h1per_wk = cum_h1per_6wk, h3per_wk = cum_h3per_6wk,
         bper_wk = cum_bper_6wk, h1per_ssn = cum_h1per,
         h3per_ssn = cum_h3per, bper_ssn = cum_bper)

# Create directory to store forecasts
dir.create("Forecasts/Training/2012-2013/Subtype Historical Average",
           showWarnings = FALSE)

# Create target densities and functions
subtype_densities_1213 <- create_subtype_densities(
  train_ili_1213, vir_ssn_per = vir_ssn_per, vir_wk_per = vir_wk_per, 
  pseudo_onsets = pseudo_onsets
)

subtype_functions_1213 <- modify_depth(
  subtype_densities_1213, 3,
  function(dens) approxfun(dens$x, dens$y, rule = 2)
)

for(i in 40:72) {
  temp <- create_subtype_forecast(
    functions = subtype_functions_1213,
    virologic = virologic_1213,
    pub_week = i,
    season = "2012/2013",
    prob_no_onset = filter(prob_no_onset, season == "2012/2013")
  )
  
  temp$Value <- format(temp$Value, scientific = FALSE)
  
  j <- str_pad(ifelse(i > 52, i - 52, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("Forecasts/Training/2012-2013/Subtype Historical Average/EW", j, ".csv"))
  
  year <- ifelse(i > 52, 2013, 2012)
  write_csv(temp,
            path = paste0(kudu_path, "/EW", j, "-", year, "-Protea_Kudu.csv"))

}

# Create forecasts for 2013/2014 ------
train_ili_1314 <- ili_current %>%
  filter(year <= 2013, season != "2013/2014",
         !location %in% state.name)

virologic_1314 <- virologic_combined %>%
  filter(season == "2013/2014", !location %in% state.name) %>%
  # Create variable to keep weeks in order
  mutate(order_week = week_inorder(week, season)) %>%
  select(year, location, week, order_week, 
         h1per_wk = cum_h1per_6wk, h3per_wk = cum_h3per_6wk,
         bper_wk = cum_bper_6wk, h1per_ssn = cum_h1per,
         h3per_ssn = cum_h3per, bper_ssn = cum_bper)

# Create directory to store forecasts
dir.create("Forecasts/Training/2013-2014/Subtype Historical Average",
           showWarnings = FALSE)

# Create target densities and functions
subtype_densities_1314 <- create_subtype_densities(
  train_ili_1314, vir_ssn_per = vir_ssn_per, vir_wk_per = vir_wk_per, 
  pseudo_onsets = pseudo_onsets
)

subtype_functions_1314 <- modify_depth(
  subtype_densities_1314, 3,
  function(dens) approxfun(dens$x, dens$y, rule = 2)
)

for(i in 40:72) {
  
  temp <- create_subtype_forecast(
    functions = subtype_functions_1314,
    virologic = virologic_1314,
    pub_week = i,
    season = "2013/2014",
    prob_no_onset = filter(prob_no_onset, season == "2013/2014")
  )
  
  temp$Value <- format(temp$Value, scientific = FALSE)
  
  j <- str_pad(ifelse(i > 52, i - 52, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("Forecasts/Training/2013-2014/Subtype Historical Average/EW", j, ".csv"))
  
  year <- ifelse(i > 52, 2014, 2013)
  write_csv(temp,
            path = paste0(kudu_path, "/EW", j, "-", year, "-Protea_Kudu.csv"))
  
}

# Create forecasts for 2014/2015 ------
train_ili_1415 <- ili_current %>%
  filter(year <= 2014, season != "2014/2015",
         !location %in% state.name)

virologic_1415 <- virologic_combined %>%
  filter(season == "2014/2015", !location %in% state.name) %>%
  # Create variable to keep weeks in order
  mutate(order_week = week_inorder(week, season)) %>%
  select(year, location, week, order_week, 
         h1per_wk = cum_h1per_6wk, h3per_wk = cum_h3per_6wk,
         bper_wk = cum_bper_6wk, h1per_ssn = cum_h1per,
         h3per_ssn = cum_h3per, bper_ssn = cum_bper)

# Create directory to store forecasts
dir.create("Forecasts/Training/2014-2015/Subtype Historical Average",
           showWarnings = FALSE)

# Create target densities and functions
subtype_densities_1415 <- create_subtype_densities(
  train_ili_1415, vir_ssn_per = vir_ssn_per, vir_wk_per = vir_wk_per, 
  pseudo_onsets = pseudo_onsets
)

subtype_functions_1415 <- modify_depth(
  subtype_densities_1415, 3,
  function(dens) approxfun(dens$x, dens$y, rule = 2)
)

for(i in 40:73) {
  temp <- create_subtype_forecast(
    functions = subtype_functions_1415,
    virologic = virologic_1415,
    pub_week = i,
    season = "2014/2015",
    prob_no_onset = filter(prob_no_onset, season == "2014/2015")
  )
  
  temp$Value <- format(temp$Value, scientific = FALSE)
  
  j <- str_pad(ifelse(i > 53, i - 53, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("Forecasts/Training/2014-2015/Subtype Historical Average/EW", j, ".csv"))
  
  year <- ifelse(i > 53, 2015, 2014)
  write_csv(temp,
            path = paste0(kudu_path, "/EW", j, "-", year, "-Protea_Kudu.csv"))
  
}

# Create forecasts for 2015/2016 ------
train_ili_1516 <- ili_current %>%
  filter(year <= 2015, season != "2015/2016",
         !location %in% state.name)

virologic_1516 <- virologic_combined %>%
  filter(season == "2015/2016", !location %in% state.name) %>%
  # Create variable to keep weeks in order
  mutate(order_week = week_inorder(week, season)) %>%
  select(year, location, week, order_week, 
         h1per_wk = cum_h1per_6wk, h3per_wk = cum_h3per_6wk,
         bper_wk = cum_bper_6wk, h1per_ssn = cum_h1per,
         h3per_ssn = cum_h3per, bper_ssn = cum_bper)

# Create directory to store forecasts
dir.create("Forecasts/Training/2015-2016/Subtype Historical Average",
           showWarnings = FALSE)

# Create target densities and functions
subtype_densities_1516 <- create_subtype_densities(
  train_ili_1516, vir_ssn_per = vir_ssn_per, vir_wk_per = vir_wk_per, 
  pseudo_onsets = pseudo_onsets
)

subtype_functions_1516 <- modify_depth(
  subtype_densities_1516, 3,
  function(dens) approxfun(dens$x, dens$y, rule = 2)
)

for(i in 40:72) {
  temp <- create_subtype_forecast(
    functions = subtype_functions_1516,
    virologic = virologic_1516,
    pub_week = i,
    season = "2015/2016",
    prob_no_onset = filter(prob_no_onset, season == "2015/2016")
  )
  
  temp$Value <- format(temp$Value, scientific = FALSE)
  
  j <- str_pad(ifelse(i > 52, i - 52, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("Forecasts/Training/2015-2016/Subtype Historical Average/EW", j, ".csv"))
  
  year <- ifelse(i > 52, 2016, 2015)
  write_csv(temp,
            path = paste0(kudu_path, "/EW", j, "-", year, "-Protea_Kudu.csv"))
  
}

# Create forecasts for 2016/2017 ------
train_ili_1617 <- ili_current %>%
  filter(year <= 2016, season != "2016/2017",
         !location %in% state.name)

virologic_1617 <- virologic_combined %>%
  filter(season == "2016/2017", !location %in% state.name) %>%
  # Create variable to keep weeks in order
  mutate(order_week = week_inorder(week, season)) %>%
  select(year, location, week, order_week, 
         h1per_wk = cum_h1per_6wk, h3per_wk = cum_h3per_6wk,
         bper_wk = cum_bper_6wk, h1per_ssn = cum_h1per,
         h3per_ssn = cum_h3per, bper_ssn = cum_bper)

# Create directory to store forecasts
dir.create("Forecasts/Training/2016-2017/Subtype Historical Average",
           showWarnings = FALSE)

# Create target densities and functions
subtype_densities_1617 <- create_subtype_densities(
  ili_df = train_ili_1617,  vir_ssn_per = vir_ssn_per, vir_wk_per = vir_wk_per, 
  pseudo_onsets = pseudo_onsets
)


subtype_functions_1617 <- modify_depth(
  subtype_densities_1617, 3,
  function(dens) approxfun(dens$x, dens$y, rule = 2)
)

for(i in 40:72) {
  temp <- create_subtype_forecast(
    functions = subtype_functions_1617,
    virologic = virologic_1617,
    pub_week = i,
    season = "2016/2017",
    prob_no_onset = filter(prob_no_onset, season == "2016/2017")
  )
  
  temp$Value <- format(temp$Value, scientific = FALSE)
  
  j <- str_pad(ifelse(i > 52, i - 52, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("Forecasts/Training/2016-2017/Subtype Historical Average/EW", j, ".csv"))
  
  year <- ifelse(i > 52, 2017, 2016)
  write_csv(temp,
            path = paste0(kudu_path, "/EW", j, "-", year, "-Protea_Kudu.csv"))
  
}

# Create forecasts for 2017/2018 ------
train_ili_1718 <- ili_current %>%
  filter(year <= 2017, season != "2017/2018",
         !location %in% state.name)

virologic_1718 <- virologic_combined %>%
  filter(season == "2017/2018", !location %in% state.name) %>%
  # Create variable to keep weeks in order
  mutate(order_week = week_inorder(week, season)) %>%
  select(year, location, week, order_week, 
         h1per_wk = cum_h1per_6wk, h3per_wk = cum_h3per_6wk,
         bper_wk = cum_bper_6wk, h1per_ssn = cum_h1per,
         h3per_ssn = cum_h3per, bper_ssn = cum_bper)

# Create directory to store forecasts
dir.create("Forecasts/Training/2017-2018/Subtype Historical Average",
           showWarnings = FALSE)

# Create target densities and functions
subtype_densities_1718 <- create_subtype_densities(
  train_ili_1718,  vir_ssn_per = vir_ssn_per, vir_wk_per = vir_wk_per, 
  pseudo_onsets = pseudo_onsets
)

subtype_functions_1718 <- modify_depth(
  subtype_densities_1718, 3,
  function(dens) approxfun(dens$x, dens$y, rule = 2)
)

for(i in 40:72) {
  temp <- create_subtype_forecast(
    functions = subtype_functions_1718,
    virologic = virologic_1718,
    pub_week = i,
    season = "2017/2018",
    prob_no_onset = filter(prob_no_onset, season == "2017/2018")
  )
  
  temp$Value <- format(temp$Value, scientific = FALSE)
  
  j <- str_pad(ifelse(i > 52, i - 52, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("Forecasts/Training/2017-2018/Subtype Historical Average/EW", j, ".csv"))
  
  year <- ifelse(i > 52, 2018, 2017)
  write_csv(temp,
            path = paste0(kudu_path, "/EW", j, "-", year, "-Protea_Kudu.csv"))
  
}

# Create forecasts for 2018/2019 ------
train_ili_1819 <- ili_current %>%
  filter(year <= 2018, season != "2018/2019",
         !location %in% state.name)

virologic_1819 <- virologic_combined %>%
  filter(season == "2018/2019", !location %in% state.name) %>%
  # Create variable to keep weeks in order
  mutate(order_week = week_inorder(week, season)) %>%
  select(year, location, week, order_week, 
         h1per_wk = cum_h1per_6wk, h3per_wk = cum_h3per_6wk,
         bper_wk = cum_bper_6wk, h1per_ssn = cum_h1per,
         h3per_ssn = cum_h3per, bper_ssn = cum_bper)

# Create directory to store forecasts
dir.create("Forecasts/Training/2018-2019/Subtype Historical Average",
           showWarnings = FALSE)

# Create target densities and functions
subtype_densities_1819 <- create_subtype_densities(
  train_ili_1819,  vir_ssn_per = vir_ssn_per, vir_wk_per = vir_wk_per, 
  pseudo_onsets = pseudo_onsets
)

subtype_functions_1819 <- modify_depth(
  subtype_densities_1819, 3,
  function(dens) approxfun(dens$x, dens$y, rule = 2)
)

for(i in 40:72) {
  temp <- create_subtype_forecast(
    functions = subtype_functions_1819,
    virologic = virologic_1819,
    pub_week = i,
    season = "2018/2019",
    prob_no_onset = filter(prob_no_onset, season == "2018/2019")
  )
  
  temp$Value <- format(temp$Value, scientific = FALSE)
  
  j <- str_pad(ifelse(i > 52, i - 52, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("Forecasts/Training/2018-2019/Subtype Historical Average/EW", j, ".csv"))
  
  year <- ifelse(i > 52, 2019, 2018)
  write_csv(temp,
            path = paste0(kudu_path, "/EW", j, "-", year, "-Protea_Kudu.csv"))
  
}