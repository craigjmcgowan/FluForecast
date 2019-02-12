# Generate historical average forecast
library(tidyverse)
library(FluSight)

# Load data -------
load("Data/ili_state.Rdata")

# Load functions 
source("R/create_historical_forecast.R")
source("R/utils.R")


# Create forecasts for 2014/2015 ------
train_ili_1415 <- state_current %>%
  filter(year <= 2014, season != "2014/2015")

# Create directory to store forecasts
dir.create("State Forecasts/2014-2015/Historical Average", 
           showWarnings = FALSE)

# Create target densities and functions
historical_densities_1415 <- create_historical_densities(train_ili_1415,
                                                         challenge = "state_ili")

historical_functions_1415 <- modify_depth(
  historical_densities_1415, 2,
  function(dens) approxfun(dens$x, dens$y, rule = 2)
)

for(i in 40:73) {
  temp <- create_historical_forecast(
    functions = historical_functions_1415,
    pub_week = i,
    season = "2014/2015",
    challenge = "state_ili"
  )
  
  j <- str_pad(ifelse(i > 53, i - 53, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("State Forecasts/2014-2015/Historical Average/EW", j, ".csv"))
  
}

# Create forecasts for 2015/2016 ------
train_ili_1516 <- state_current %>%
  filter(year <= 2015, season != "2015/2016")

# Create directory to store forecasts
dir.create("State Forecasts/2015-2016/Historical Average",
           showWarnings = FALSE)

# Create target densities and functions
historical_densities_1516 <- create_historical_densities(train_ili_1516,
                                                         challenge = "state_ili")

historical_functions_1516 <- modify_depth(
  historical_densities_1516, 2,
  function(dens) approxfun(dens$x, dens$y, rule = 2)
)

for(i in 40:72) {
  temp <- create_historical_forecast(
    functions = historical_functions_1516,
    pub_week = i,
    season = "2015/2016",
    challenge = "state_ili"
  )
  
  j <- str_pad(ifelse(i > 52, i - 52, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("State Forecasts/2015-2016/Historical Average/EW", j, ".csv"))
  
}

# Create forecasts for 2016/2017 ------
train_ili_1617 <- state_current %>%
  filter(year <= 2016, season != "2016/2017")

# Create directory to store forecasts
dir.create("State Forecasts/2016-2017/Historical Average",
           showWarnings = FALSE)

# Create target densities and functions
historical_densities_1617 <- create_historical_densities(train_ili_1617,
                                                         challenge = "state_ili")

historical_functions_1617 <- modify_depth(
  historical_densities_1617, 2,
  function(dens) approxfun(dens$x, dens$y, rule = 2)
)

for(i in 40:72) {
  temp <- create_historical_forecast(
    functions = historical_functions_1617,
    pub_week = i,
    season = "2016/2017",
    challenge = "state_ili"
  )
  
  j <- str_pad(ifelse(i > 52, i - 52, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("State Forecasts/2016-2017/Historical Average/EW", j, ".csv"))
  
}

# Create forecasts for 2017/2018 ------
train_ili_1718 <- state_current %>%
  filter(year <= 2017, season != "2017/2018")

# Create directory to store forecasts
dir.create("State Forecasts/2017-2018/Historical Average",
           showWarnings = FALSE)

# Create target densities and functions
historical_densities_1718 <- create_historical_densities(train_ili_1718,
                                                         challenge = "state_ili")

historical_functions_1718 <- modify_depth(
  historical_densities_1718, 2,
  function(dens) approxfun(dens$x, dens$y, rule = 2)
)

for(i in 40:72) {
  temp <- create_historical_forecast(
    functions = historical_functions_1718,
    pub_week = i,
    season = "2017/2018",
    challenge = "state_ili"
  )
  
  j <- str_pad(ifelse(i > 52, i - 52, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("State Forecasts/2017-2018/Historical Average/EW", j, ".csv"))
  
}

# Create forecasts for 2018/2019 ------
train_ili_1819 <- state_current %>%
  filter(year <= 2018, season != "2018/2019")

# Create directory to store forecasts
dir.create("State Forecasts/2018-2019/Historical Average",
           showWarnings = FALSE)

# Create target densities and functions
historical_densities_1819 <- create_historical_densities(train_ili_1819,
                                                         challenge = "state_ili")

historical_functions_1819 <- modify_depth(
  historical_densities_1819, 2,
  function(dens) approxfun(dens$x, dens$y, rule = 2)
)

for(i in 40:72) {
  temp <- create_historical_forecast(
    functions = historical_functions_1819,
    pub_week = i,
    season = "2018/2019",
    challenge = "state_ili"
  )
  
  j <- str_pad(ifelse(i > 52, i - 52, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("State Forecasts/2018-2019/Historical Average/EW", j, ".csv"))
  
}
