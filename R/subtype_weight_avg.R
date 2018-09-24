# Generate subtype-weighted historical average forecast
library(tidyverse)
library(FluSight)
library(zoo)

# Load data -------
load("Data/ili.Rdata")
load("Data/virologic.Rdata")
load("Data/pseudo_onsets_2003_2006.Rdata")

# Load functions 
source("R/create_subtype_forecast.R")
source("R/utils.R")

# Calculate cumulative virologic status across seasons -----
vir_ssn_per <- virologic_combined %>%
  group_by(season, location) %>%
  summarize(h1sum = sum(a_2009_h1n1, na.rm = T) + sum(a_h1, na.rm = T),
            h3sum = sum(a_h3, na.rm = T),
            h1per = h1sum / (h1sum + h3sum),
            h3per = h3sum / (h1sum + h3sum)) %>%
  select(-h1sum, -h3sum)

# Create data.frame of probability of no onset based on prior years -----
onsets <- ili_current %>%
  filter(year >= 2007, !season %in% c("2006/2007", "2009/2010"),
         (week >= 40 | week <= 20)) %>%
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
         h3_prob_no_onset = ifelse(isna, NA, cummean(h3per))) %>%
  group_by(location) %>%
  mutate(h1_prob_no_onset = na.locf(h1_prob_no_onset, na.rm = FALSE),
         h3_prob_no_onset = na.locf(h3_prob_no_onset, na.rm = FALSE)) %>%
  mutate(h1_prob_no_onset = ifelse(is.na(h1_prob_no_onset), 0, 
                                   h1_prob_no_onset * prob_no_onset),
         h3_prob_no_onset = ifelse(is.na(h3_prob_no_onset), 0, 
                                   h3_prob_no_onset * prob_no_onset)) %>%
  ungroup() %>%
  select(location, season, prob_no_onset, h1_prob_no_onset, h3_prob_no_onset)

# Create forecasts for 2010/2011 ------
train_ili_1011 <- ili_current %>%
  filter(year <= 2010, season != "2010/2011")

virologic_1011 <- virologic_combined %>%
  filter(season == "2010/2011") %>%
  group_by(location) %>%
  mutate(h1sum = cumsum(a_2009_h1n1) + cumsum(a_h1),
         h3sum = cumsum(a_h3),
         h1per = h1sum / (h1sum + h3sum),
         h3per = h3sum / (h1sum + h3sum)) %>%
  # If no samples reported, make each subtype 50%
  mutate(h1per = ifelse(is.na(h1per), 0.5, h1per),
         h3per = ifelse(is.na(h3per), 0.5, h3per)) %>%
  # Create variable to keep weeks in order
  mutate(order_week = week_inorder(week, season)) %>%
  select(year,location, week, h1sum, h1per, h3sum, h3per, order_week)

# Create directory to store forecasts
dir.create("Forecasts/2010-2011/Subtype Historical Average/",
           showWarnings = FALSE)

# Create target densities and functions
subtype_densities_1011 <- create_subtype_densities(
  train_ili_1011, vir_ssn_per, pseudo_onsets
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
  
  j <- str_pad(ifelse(i > 52, i - 52, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("Forecasts/2010-2011/Subtype Historical Average/EW", j, ".csv"))
  
}

# Create forecasts for 2011/2012 ------
train_ili_1112 <- ili_current %>%
  filter(year <= 2011, season != "2011/2012")

virologic_1112 <- virologic_combined %>%
  filter(season == "2011/2012") %>%
  group_by(location) %>%
  mutate(h1sum = cumsum(a_2009_h1n1) + cumsum(a_h1),
         h3sum = cumsum(a_h3),
         h1per = h1sum / (h1sum + h3sum),
         h3per = h3sum / (h1sum + h3sum)) %>%
  # If no samples reported, make each subtype 50%
  mutate(h1per = ifelse(is.na(h1per), 0.5, h1per),
         h3per = ifelse(is.na(h3per), 0.5, h3per)) %>%
  # Create variable to keep weeks in order
  mutate(order_week = week_inorder(week, season)) %>%
  select(year,location, week, h1sum, h1per, h3sum, h3per, order_week)

# Create directory to store forecasts
dir.create("Forecasts/2011-2012/Subtype Historical Average/",
           showWarnings = FALSE)

# Create target densities and functions
subtype_densities_1112 <- create_subtype_densities(
  train_ili_1112, vir_ssn_per, pseudo_onsets
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
  
  j <- str_pad(ifelse(i > 52, i - 52, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("Forecasts/2011-2012/Subtype Historical Average/EW", j, ".csv"))
  
}

# Create forecasts for 2012/2013 ------
train_ili_1213 <- ili_current %>%
  filter(year <= 2012, season != "2012/2013")

virologic_1213 <- virologic_combined %>%
  filter(season == "2012/2013") %>%
  group_by(location) %>%
  mutate(h1sum = cumsum(a_2009_h1n1) + cumsum(a_h1),
         h3sum = cumsum(a_h3),
         h1per = h1sum / (h1sum + h3sum),
         h3per = h3sum / (h1sum + h3sum)) %>%
  # If no samples reported, make each subtype 50%
  mutate(h1per = ifelse(is.na(h1per), 0.5, h1per),
         h3per = ifelse(is.na(h3per), 0.5, h3per)) %>%
  # Create variable to keep weeks in order
  mutate(order_week = week_inorder(week, season)) %>%
  select(year,location, week, h1sum, h1per, h3sum, h3per, order_week)

# Create directory to store forecasts
dir.create("Forecasts/2012-2013/Subtype Historical Average",
           showWarnings = FALSE)

# Create target densities and functions
subtype_densities_1213 <- create_subtype_densities(
  train_ili_1213, vir_ssn_per, pseudo_onsets
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
  
  j <- str_pad(ifelse(i > 52, i - 52, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("Forecasts/2012-2013/Subtype Historical Average/EW", j, ".csv"))

}

# Create forecasts for 2013/2014 ------
train_ili_1314 <- ili_current %>%
  filter(year <= 2014, season != "2014/2015")

virologic_1314 <- virologic_combined %>%
  filter(season == "2013/2014") %>%
  group_by(location) %>%
  mutate(h1sum = cumsum(a_2009_h1n1) + cumsum(a_h1),
         h3sum = cumsum(a_h3),
         h1per = h1sum / (h1sum + h3sum),
         h3per = h3sum / (h1sum + h3sum)) %>%
  # If no samples reported, make each subtype 50%
  mutate(h1per = ifelse(is.na(h1per), 0.5, h1per),
         h3per = ifelse(is.na(h3per), 0.5, h3per)) %>%
  # Create variable to keep weeks in order
  mutate(order_week = week_inorder(week, season)) %>%
  select(year,location, week, h1sum, h1per, h3sum, h3per, order_week)

# Create directory to store forecasts
dir.create("Forecasts/2013-2014/Subtype Historical Average",
           showWarnings = FALSE)

# Create target densities and functions
subtype_densities_1314 <- create_subtype_densities(
  train_ili_1314, vir_ssn_per, pseudo_onsets
)
subtype_functions_1314 <- modify_depth(
  subtype_densities_1314, 3,
  function(dens) approxfun(dens$x, dens$y, rule = 2)
)

for(i in 40:72) {
  j <- str_pad(ifelse(i > 52, i - 52, i), 2, pad = "0")
  
  write_csv(create_subtype_forecast(
      functions = subtype_functions_1314,
      virologic = virologic_1314,
      pub_week = i,
      season = "2013/2014",
      prob_no_onset = filter(prob_no_onset, season == "2013/2014")
    ), path = paste0("Forecasts/2013-2014/Subtype Historical Average/EW", j, ".csv"))
  
}

# Create forecasts for 2014/2015 ------
train_ili_1415 <- ili_current %>%
  filter(year <= 2015, season != "2015/2016")

virologic_1415 <- virologic_combined %>%
  filter(season == "2014/2015") %>%
  group_by(location) %>%
  mutate(h1sum = cumsum(a_2009_h1n1) + cumsum(a_h1),
         h3sum = cumsum(a_h3),
         h1per = h1sum / (h1sum + h3sum),
         h3per = h3sum / (h1sum + h3sum)) %>%
  # If no samples reported, make each subtype 50%
  mutate(h1per = ifelse(is.na(h1per), 0.5, h1per),
         h3per = ifelse(is.na(h3per), 0.5, h3per)) %>%
  # Create variable to keep weeks in order
  mutate(order_week = week_inorder(week, season)) %>%
  select(year,location, week, h1sum, h1per, h3sum, h3per, order_week)

# Create directory to store forecasts
dir.create("Forecasts/2014-2015/Subtype Historical Average",
           showWarnings = FALSE)

# Create target densities and functions
subtype_densities_1415 <- create_subtype_densities(
  train_ili_1415, vir_ssn_per, pseudo_onsets
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
  
  j <- str_pad(ifelse(i > 53, i - 53, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("Forecasts/2014-2015/Subtype Historical Average/EW", j, ".csv"))
  
}

# Create forecasts for 2015/2016 ------
train_ili_1516 <- ili_current %>%
  filter(year <= 2015, season != "2015/2016")

virologic_1516 <- virologic_combined %>%
  filter(season == "2015/2016") %>%
  group_by(location) %>%
  mutate(h1sum = cumsum(a_2009_h1n1) + cumsum(a_h1),
         h3sum = cumsum(a_h3),
         h1per = h1sum / (h1sum + h3sum),
         h3per = h3sum / (h1sum + h3sum)) %>%
  # If no samples reported, make each subtype 50%
  mutate(h1per = ifelse(is.na(h1per), 0.5, h1per),
         h3per = ifelse(is.na(h3per), 0.5, h3per)) %>%
  # Create variable to keep weeks in order
  mutate(order_week = week_inorder(week, season)) %>%
  select(year,location, week, h1sum, h1per, h3sum, h3per, order_week)

# Create directory to store forecasts
dir.create("Forecasts/2015-2016/Subtype Historical Average",
           showWarnings = FALSE)

# Create target densities and functions
subtype_densities_1516 <- create_subtype_densities(
  train_ili_1516, vir_ssn_per, pseudo_onsets
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
  
  j <- str_pad(ifelse(i > 52, i - 52, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("Forecasts/2015-2016/Subtype Historical Average/EW", j, ".csv"))
  
}

# Create forecasts for 2016/2017 ------
train_ili_1617 <- ili_current %>%
  filter(year <= 2016, season != "2016/2017")

virologic_1617 <- virologic_combined %>%
  filter(season == "2016/2017") %>%
  group_by(location) %>%
  mutate(h1sum = cumsum(a_2009_h1n1) + cumsum(a_h1),
         h3sum = cumsum(a_h3),
         h1per = h1sum / (h1sum + h3sum),
         h3per = h3sum / (h1sum + h3sum)) %>%
  # If no samples reported, make each subtype 50%
  mutate(h1per = ifelse(is.na(h1per), 0.5, h1per),
         h3per = ifelse(is.na(h3per), 0.5, h3per)) %>%
  # Create variable to keep weeks in order
  mutate(order_week = week_inorder(week, season)) %>%
  select(year,location, week, h1sum, h1per, h3sum, h3per, order_week)

# Create directory to store forecasts
dir.create("Forecasts/2016-2017/Subtype Historical Average",
           showWarnings = FALSE)

# Create target densities and functions
subtype_densities_1617 <- create_subtype_densities(
  ili_df = train_ili_1617, vir_ssn_per = vir_ssn_per, pseudo_onsets
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
  
  j <- str_pad(ifelse(i > 52, i - 52, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("Forecasts/2016-2017/Subtype Historical Average/EW", j, ".csv"))
  
}

# Create forecasts for 2017/2018 ------
train_ili_1718 <- ili_current %>%
  filter(year <= 2018, season != "2017/2018")

virologic_1718 <- virologic_combined %>%
  filter(season == "2017/2018") %>%
  group_by(location) %>%
  mutate(h1sum = cumsum(a_2009_h1n1) + cumsum(a_h1),
         h3sum = cumsum(a_h3),
         h1per = h1sum / (h1sum + h3sum),
         h3per = h3sum / (h1sum + h3sum)) %>%
  # If no samples reported, make each subtype 50%
  mutate(h1per = ifelse(is.na(h1per), 0.5, h1per),
         h3per = ifelse(is.na(h3per), 0.5, h3per)) %>%
  # Create variable to keep weeks in order
  mutate(order_week = week_inorder(week, season)) %>%
  select(year,location, week, h1sum, h1per, h3sum, h3per, order_week)

# Create directory to store forecasts
dir.create("Forecasts/2017-2018/Subtype Historical Average",
           showWarnings = FALSE)

# Create target densities and functions
subtype_densities_1718 <- create_subtype_densities(
  train_ili_1718, vir_ssn_per, pseudo_onsets
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
  
  j <- str_pad(ifelse(i > 52, i - 52, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("Forecasts/2017-2018/Subtype Historical Average/EW", j, ".csv"))
  
}