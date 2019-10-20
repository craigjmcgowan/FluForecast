# Generate historical average forecast
library(tidyverse)
library(FluSight)

# Load data -------
ili_current <- readRDS("Data/ili_current.RDS")
pseudo_onsets <- readRDS("Data/pseudo_onsets.RDS")

# Load functions 
source("R/create_historical_forecast.R")
source("R/utils.R")

# Create data.frame of probability of no onset based on prior years -----
onsets <- ili_current %>%
  filter(year >= 2007, season != "2006/2007",
         (week >= 40 | week <= 20),
         !location %in% state.name) %>%
  group_by(location, season) %>%
  do(create_onset(., region = .$location[1], 
                  year = as.numeric(substr(.$season[1], 1, 4)))) %>%
  ungroup() 

# Weight probability of no onset by virus % in seasons w/o onset
prob_no_onset <- onsets %>%
  group_by(location) %>%
  arrange(season, .by_group = T) %>%
  mutate(prob_no_onset = lag(cumsum(bin_start_incl == "none") / row_number(), default = 0))


# Create forecasts for 2010/2011 ------
train_ili_1011 <- ili_current %>%
  filter(year <= 2010, season != "2010/2011",
         !location %in% state.name)

# Create directory to store forecasts
dir.create("Forecasts/2010-2011/Historical Average",
           showWarnings = FALSE)

# Create target densities and functions
historical_densities_1011 <- create_historical_densities(train_ili_1011,
                                                         pseudo_onsets)

historical_functions_1011 <- modify_depth(
  historical_densities_1011, 2,
  function(dens) approxfun(dens$x, dens$y, rule = 2)
)

for(i in 40:72) {
  temp <- create_historical_forecast(
    functions = historical_functions_1011,
    pub_week = i,
    season = "2010/2011",
    prob_no_onset = filter(prob_no_onset, season == "2010/2011")
  )
  
  j <- str_pad(ifelse(i > 52, i - 52, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("Forecasts/2010-2011/Historical Average/EW", j, ".csv"))
  
}

# Create forecasts for 2011/2012 ------
train_ili_1112 <- ili_current %>%
  filter(year <= 2011, season != "2011/2012",
         !location %in% state.name)

# Create directory to store forecasts
dir.create("Forecasts/2011-2012/Historical Average",
           showWarnings = FALSE)

# Create target densities and functions
historical_densities_1112 <- create_historical_densities(train_ili_1112,
                                                         pseudo_onsets)

historical_functions_1112 <- modify_depth(
  historical_densities_1112, 2,
  function(dens) approxfun(dens$x, dens$y, rule = 2)
)

for(i in 40:72) {
  temp <- create_historical_forecast(
    functions = historical_functions_1112,
    pub_week = i,
    season = "2011/2012",
    prob_no_onset = filter(prob_no_onset, season == "2011/2012")
  )
  
  j <- str_pad(ifelse(i > 52, i - 52, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("Forecasts/2011-2012/Historical Average/EW", j, ".csv"))
  
}


# Create forecasts for 2012/2013 ------
train_ili_1213 <- ili_current %>%
  filter(year <= 2012, season != "2012/2013",
         !location %in% state.name)

# Create directory to store forecasts
dir.create("Forecasts/2012-2013/Historical Average",
           showWarnings = FALSE)

# Create target densities and functions
historical_densities_1213 <- create_historical_densities(train_ili_1213,
                                                         pseudo_onsets)

historical_functions_1213 <- modify_depth(
  historical_densities_1213, 2,
  function(dens) approxfun(dens$x, dens$y, rule = 2)
)

for(i in 40:72) {
  temp <- create_historical_forecast(
    functions = historical_functions_1213,
    pub_week = i,
    season = "2012/2013",
    prob_no_onset = filter(prob_no_onset, season == "2012/2013")
  )
  
  j <- str_pad(ifelse(i > 52, i - 52, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("Forecasts/2012-2013/Historical Average/EW", j, ".csv"))

}

# Create forecasts for 2013/2014 ------
train_ili_1314 <- ili_current %>%
  filter(year <= 2013, season != "2013/2014",
         !location %in% state.name)

# Create directory to store forecasts
dir.create("Forecasts/2013-2014/Historical Average/",
           showWarnings = FALSE)

# Create target densities and functions
historical_densities_1314 <- create_historical_densities(train_ili_1314,
                                                         pseudo_onsets)

historical_functions_1314 <- modify_depth(
  historical_densities_1314, 2,
  function(dens) approxfun(dens$x, dens$y, rule = 2)
)

for(i in 40:72) {
  temp <- create_historical_forecast(
    functions = historical_functions_1314,
    pub_week = i,
    season = "2013/2014",
    prob_no_onset = filter(prob_no_onset, season == "2013/2014")
  )
  
  j <- str_pad(ifelse(i > 52, i - 52, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("Forecasts/2013-2014/Historical Average/EW", j, ".csv"))
  
}

# Create forecasts for 2014/2015 ------
train_ili_1415 <- ili_current %>%
  filter(year <= 2014, season != "2014/2015",
         !location %in% state.name)

# Create directory to store forecasts
dir.create("Forecasts/2014-2015/Historical Average",
           showWarnings = FALSE)

# Create target densities and functions
historical_densities_1415 <- create_historical_densities(train_ili_1415,
                                                         pseudo_onsets)

historical_functions_1415 <- modify_depth(
  historical_densities_1415, 2,
  function(dens) approxfun(dens$x, dens$y, rule = 2)
)

for(i in 40:73) {
  temp <- create_historical_forecast(
    functions = historical_functions_1415,
    pub_week = i,
    season = "2014/2015",
    prob_no_onset = filter(prob_no_onset, season == "2014/2015")
  )
  
  j <- str_pad(ifelse(i > 53, i - 53, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("Forecasts/2014-2015/Historical Average/EW", j, ".csv"))
  
}

# State forecasts
state_train_ili_1415 <- ili_current %>%
  filter(year <= 2014, season != "2014/2015",
         location %in% state.name)

# Create directory to store forecasts
dir.create("State Forecasts/2014-2015/Historical Average", 
           showWarnings = FALSE)

# Create target densities and functions
state_historical_densities_1415 <- create_historical_densities(state_train_ili_1415,
                                                               challenge = "state_ili")

state_historical_functions_1415 <- modify_depth(
  state_historical_densities_1415, 2,
  function(dens) approxfun(dens$x, dens$y, rule = 2)
)

for(i in 40:73) {
  temp <- create_historical_forecast(
    functions = state_historical_functions_1415,
    pub_week = i,
    season = "2014/2015",
    challenge = "state_ili"
  )
  
  j <- str_pad(ifelse(i > 53, i - 53, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("State Forecasts/2014-2015/Historical Average/EW", j, ".csv"))
  
}


# Create forecasts for 2015/2016 ------
train_ili_1516 <- ili_current %>%
  filter(year <= 2015, season != "2015/2016",
         !location %in% state.name)

# Create directory to store forecasts
dir.create("Forecasts/2015-2016/Historical Average",
           showWarnings = FALSE)

# Create target densities and functions
historical_densities_1516 <- create_historical_densities(train_ili_1516,
                                                         pseudo_onsets)

historical_functions_1516 <- modify_depth(
  historical_densities_1516, 2,
  function(dens) approxfun(dens$x, dens$y, rule = 2)
)

for(i in 40:72) {
  temp <- create_historical_forecast(
    functions = historical_functions_1516,
    pub_week = i,
    season = "2015/2016",
    prob_no_onset = filter(prob_no_onset, season == "2015/2016")
  )
  
  j <- str_pad(ifelse(i > 52, i - 52, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("Forecasts/2015-2016/Historical Average/EW", j, ".csv"))
  
}

# State forecasts
state_train_ili_1516 <- ili_current %>%
  filter(year <= 2015, season != "2015/2016",
         location %in% state.name)

# Create directory to store forecasts
dir.create("State Forecasts/2015-2016/Historical Average",
           showWarnings = FALSE)

# Create target densities and functions
state_historical_densities_1516 <- create_historical_densities(train_ili_1516,
                                                               challenge = "state_ili")

state_historical_functions_1516 <- modify_depth(
  state_historical_densities_1516, 2,
  function(dens) approxfun(dens$x, dens$y, rule = 2)
)

for(i in 40:72) {
  temp <- create_historical_forecast(
    functions = state_historical_functions_1516,
    pub_week = i,
    season = "2015/2016",
    challenge = "state_ili"
  )
  
  j <- str_pad(ifelse(i > 52, i - 52, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("State Forecasts/2015-2016/Historical Average/EW", j, ".csv"))
  
}

# Create forecasts for 2016/2017 ------
train_ili_1617 <- ili_current %>%
  filter(year <= 2016, season != "2016/2017",
         !location %in% state.name)

# Create directory to store forecasts
dir.create("Forecasts/2016-2017/Historical Average",
           showWarnings = FALSE)

# Create target densities and functions
historical_densities_1617 <- create_historical_densities(train_ili_1617,
                                                         pseudo_onsets)

historical_functions_1617 <- modify_depth(
  historical_densities_1617, 2,
  function(dens) approxfun(dens$x, dens$y, rule = 2)
)

for(i in 40:72) {
  temp <- create_historical_forecast(
    functions = historical_functions_1617,
    pub_week = i,
    season = "2016/2017",
    prob_no_onset = filter(prob_no_onset, season == "2016/2017")
  )
  
  j <- str_pad(ifelse(i > 52, i - 52, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("Forecasts/2016-2017/Historical Average/EW", j, ".csv"))
  
}

# State forecasts
state_train_ili_1617 <- ili_current %>%
  filter(year <= 2016, season != "2016/2017",
         location %in% state.name)

# Create directory to store forecasts
dir.create("State Forecasts/2016-2017/Historical Average",
           showWarnings = FALSE)

# Create target densities and functions
state_historical_densities_1617 <- create_historical_densities(state_train_ili_1617,
                                                               challenge = "state_ili")

state_historical_functions_1617 <- modify_depth(
  state_historical_densities_1617, 2,
  function(dens) approxfun(dens$x, dens$y, rule = 2)
)

for(i in 40:72) {
  temp <- create_historical_forecast(
    functions = state_historical_functions_1617,
    pub_week = i,
    season = "2016/2017",
    challenge = "state_ili"
  )
  
  j <- str_pad(ifelse(i > 52, i - 52, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("State Forecasts/2016-2017/Historical Average/EW", j, ".csv"))
  
}

# Create forecasts for 2017/2018 ------
train_ili_1718 <- ili_current %>%
  filter(year <= 2017, season != "2017/2018",
         !location %in% state.name)

# Create directory to store forecasts
dir.create("Forecasts/2017-2018/Historical Average",
           showWarnings = FALSE)

# Create target densities and functions
historical_densities_1718 <- create_historical_densities(train_ili_1718,
                                                         pseudo_onsets)

historical_functions_1718 <- modify_depth(
  historical_densities_1718, 2,
  function(dens) approxfun(dens$x, dens$y, rule = 2)
)

for(i in 40:72) {
  temp <- create_historical_forecast(
    functions = historical_functions_1718,
    pub_week = i,
    season = "2017/2018",
    prob_no_onset = filter(prob_no_onset, season == "2017/2018")
  )
  
  j <- str_pad(ifelse(i > 52, i - 52, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("Forecasts/2017-2018/Historical Average/EW", j, ".csv"))
  
}

# State forecasts
state_train_ili_1718 <- ili_current %>%
  filter(year <= 2017, season != "2017/2018",
         location %in% state.name, location != "Louisiana") # Remove LA b/c only one year of training data

# Create directory to store forecasts
dir.create("State Forecasts/2017-2018/Historical Average",
           showWarnings = FALSE)

# Create target densities and functions
state_historical_densities_1718 <- create_historical_densities(state_train_ili_1718,
                                                               challenge = "state_ili")

state_historical_functions_1718 <- modify_depth(
  state_historical_densities_1718, 2,
  function(dens) approxfun(dens$x, dens$y, rule = 2)
)

for(i in 40:72) {
  temp <- create_historical_forecast(
    functions = state_historical_functions_1718,
    pub_week = i,
    season = "2017/2018",
    challenge = "state_ili"
  )
  
  j <- str_pad(ifelse(i > 52, i - 52, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("State Forecasts/2017-2018/Historical Average/EW", j, ".csv"))
  
}

# Create forecasts for 2018/2019 ------
train_ili_1819 <- ili_current %>%
  filter(year <= 2018, season != "2018/2019",
         !location %in% state.name)

# Create directory to store forecasts
dir.create("Forecasts/2018-2019/Historical Average",
           showWarnings = FALSE)

# Create target densities and functions
historical_densities_1819 <- create_historical_densities(train_ili_1819,
                                                         pseudo_onsets)

historical_functions_1819 <- modify_depth(
  historical_densities_1819, 2,
  function(dens) approxfun(dens$x, dens$y, rule = 2)
)

for(i in 40:72) {
  temp <- create_historical_forecast(
    functions = historical_functions_1819,
    pub_week = i,
    season = "2018/2019",
    prob_no_onset = filter(prob_no_onset, season == "2018/2019")
  )
  
  j <- str_pad(ifelse(i > 52, i - 52, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("Forecasts/2018-2019/Historical Average/EW", j, ".csv"))
  
}

# State forecasts
state_train_ili_1819 <- ili_current %>%
  filter(year <= 2018, season != "2018/2019",
         location %in% state.name)

# Create directory to store forecasts
dir.create("State Forecasts/2018-2019/Historical Average",
           showWarnings = FALSE)

# Create target densities and functions
state_historical_densities_1819 <- create_historical_densities(state_train_ili_1819,
                                                               challenge = "state_ili")

state_historical_functions_1819 <- modify_depth(
  state_historical_densities_1819, 2,
  function(dens) approxfun(dens$x, dens$y, rule = 2)
)

for(i in 40:72) {
  temp <- create_historical_forecast(
    functions = state_historical_functions_1819,
    pub_week = i,
    season = "2018/2019",
    challenge = "state_ili"
  )
  
  j <- str_pad(ifelse(i > 52, i - 52, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("State Forecasts/2018-2019/Historical Average/EW", j, ".csv"))
  
}


# Create forecasts for 2019/2020 ------
train_ili_1920<- ili_current %>%
  filter(year <= 2019, season != "2019/2020",
         !location %in% state.name)

# Create directory to store forecasts
dir.create("Forecasts/Training/2019-2020/Historical Average",
           showWarnings = FALSE)

# Create target densities and functions
historical_densities_1920 <- create_historical_densities(train_ili_1920,
                                                         pseudo_onsets)

historical_functions_1920 <- modify_depth(
  historical_densities_1920, 2,
  function(dens) approxfun(dens$x, dens$y, rule = 2)
)

for(i in 40:72) {
  temp <- create_historical_forecast(
    functions = historical_functions_1920,
    pub_week = i,
    season = "2019/2020",
    prob_no_onset = filter(prob_no_onset, season == "2019/2020")
  )
  
  j <- str_pad(ifelse(i > 52, i - 52, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("Forecasts/Live/2019-2020/Historical Average/EW", j, ".csv"))
  
}

# State forecasts
state_train_ili_1920 <- ili_current %>%
  filter(year <= 2019, season != "2019/2020",
         location %in% state.name)

# Create directory to store forecasts
dir.create("State Forecasts/2019-2020/Historical Average",
           showWarnings = FALSE)

# Create target densities and functions
state_historical_densities_1920 <- create_historical_densities(state_train_ili_1920,
                                                               challenge = "state_ili")

state_historical_functions_1920 <- modify_depth(
  state_historical_densities_1920, 2,
  function(dens) approxfun(dens$x, dens$y, rule = 2)
)

for(i in 40:72) {
  temp <- create_historical_forecast(
    functions = state_historical_functions_1920,
    pub_week = i,
    season = "2019/2020",
    challenge = "state_ili"
  )
  
  j <- str_pad(ifelse(i > 52, i - 52, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("State Forecasts/2019-2020/Historical Average/EW", j, ".csv"))
  
}