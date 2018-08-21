# Generate subtype-weighted historical average forecast
library(tidyverse)
library(FluSight)

# Load data -------
load("Data/ili.Rdata")

# Load functions 
source("R/create_historical_forecast.R")
source("R/week_order_functions.R")

# Create data.frame of probability of no onset based on prior years -----
onsets <- ili_current %>%
  filter(year >= 2007, season != "2006/2007",
         (week >= 40 | week <= 20)) %>%
  group_by(location, season) %>%
  do(create_onset(., region = .$location[1], 
                  year = as.numeric(substr(.$season[1], 1, 4)))) %>%
  ungroup() 

# Weight probability of no onset by virus % in seasons w/o onset
prob_no_onset <- onsets %>%
  group_by(location) %>%
  summarize(prob_no_onset = sum(bin_start_incl == "none") / n())


# Create forecasts for 2010/2011 ------
## Errors in generating densities for season onset - training data too sparse

# Create forecasts for 2011/2012 ------
## Errors in generating densities for season onset - training data too sparse

# Create forecasts for 2012/2013 ------
## Errors in generating densities for season onset - training data too sparse

train_ili_1213 <- ili_current %>%
  filter(year <= 2012, season != "2012/2013")

# Create directory to store forecasts
dir.create("Forecasts/2012-2013/Historical Average",
           showWarnings = FALSE)

# Create target densities and functions
historical_densities_1213 <- create_historical_densities(train_ili_1213)

historical_functions_1213 <- modify_depth(
  historical_densities_1213, 2,
  function(dens) approxfun(dens$x, dens$y, rule = 2)
)

for(i in 40:72) {
  temp <- create_historical_forecast(
    functions = historical_functions_1213,
    pub_week = i,
    season = "2012/2013",
    prob_no_onset = prob_no_onset
  )
  
  j <- str_pad(ifelse(i > 52, i - 52, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("Forecasts/2012-2013/Historical Average/EW", j, ".csv"))

}

# Create forecasts for 2013/2014 ------
train_ili_1314 <- ili_current %>%
  filter(year <= 2013, season != "2013/2014")

# Create directory to store forecasts
dir.create("Forecasts/2013-2014/Historical Average/",
           showWarnings = FALSE)

# Create target densities and functions
historical_densities_1314 <- create_historical_densities(train_ili_1314)

historical_functions_1314 <- modify_depth(
  historical_densities_1314, 2,
  function(dens) approxfun(dens$x, dens$y, rule = 2)
)

for(i in 40:72) {
  temp <- create_historical_forecast(
    functions = historical_functions_1314,
    pub_week = i,
    season = "2013/2014",
    prob_no_onset = prob_no_onset
  )
  
  j <- str_pad(ifelse(i > 52, i - 52, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("Forecasts/2013-2014/Historical Average/EW", j, ".csv"))
  
}

# Create forecasts for 2014/2015 ------
train_ili_1415 <- ili_current %>%
  filter(year <= 2014, season != "2014/2015")

# Create directory to store forecasts
dir.create("Forecasts/2014-2015/Historical Average",
           showWarnings = FALSE)

# Create target densities and functions
historical_densities_1415 <- create_historical_densities(train_ili_1415)

historical_functions_1415 <- modify_depth(
  historical_densities_1415, 2,
  function(dens) approxfun(dens$x, dens$y, rule = 2)
)

for(i in 40:73) {
  temp <- create_historical_forecast(
    functions = historical_functions_1415,
    pub_week = i,
    season = "2014/2015",
    prob_no_onset = prob_no_onset
  )
  
  j <- str_pad(ifelse(i > 52, i - 52, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("Forecasts/2014-2015/Historical Average/EW", j, ".csv"))
  
}

# Create forecasts for 2015/2016 ------
train_ili_1516 <- ili_current %>%
  filter(year <= 2015, season != "2015/2016")

# Create directory to store forecasts
dir.create("Forecasts/Historical Average/2015-2016",
           showWarnings = FALSE)

# Create target densities and functions
historical_densities_1516 <- create_historical_densities(train_ili_1516)

historical_functions_1516 <- modify_depth(
  historical_densities_1516, 2,
  function(dens) approxfun(dens$x, dens$y, rule = 2)
)

for(i in 40:72) {
  temp <- create_historical_forecast(
    functions = historical_functions_1516,
    pub_week = i,
    season = "2015/2016",
    prob_no_onset = prob_no_onset
  )
  
  j <- str_pad(ifelse(i > 52, i - 52, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("Forecasts/Historical Average/2015-2016/EW", j, ".csv"))
  
}

# Create forecasts for 2016/2017 ------
train_ili_1617 <- ili_current %>%
  filter(year <= 2017, season != "2017/2018")

# Create directory to store forecasts
dir.create("Forecasts/Historical Average/2016-2017",
           showWarnings = FALSE)

# Create target densities and functions
historical_densities_1617 <- create_historical_densities(train_ili_1617)

historical_functions_1617 <- modify_depth(
  historical_densities_1617, 1,
  function(dens) approxfun(dens$x, dens$y, rule = 2)
)

for(i in 40:72) {
  temp <- create_historical_forecast(
    functions = historical_functions_1617,
    pub_week = i,
    season = "2016/2017",
    prob_no_onset = prob_no_onset
  )
  
  j <- str_pad(ifelse(i > 52, i - 52, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("Forecasts/Historical Average/2016-2017/EW", j, ".csv"))
  
}

# Create forecasts for 2017/2018 ------
train_ili_1718 <- ili_current %>%
  filter(year <= 2018, season != "2018/2019")


# Create directory to store forecasts
dir.create("Forecasts/Historical Average/2017-2018",
           showWarnings = FALSE)

# Create target densities and functions
historical_densities_1718 <- create_historical_densities(train_ili_1718)

historical_functions_1718 <- modify_depth(
  historical_densities_1718, 2,
  function(dens) approxfun(dens$x, dens$y, rule = 2)
)

for(i in 40:72) {
  temp <- create_historical_forecast(
    functions = historical_functions_1718,
    pub_week = i,
    season = "2017/2018",
    prob_no_onset = prob_no_onset
  )
  
  j <- str_pad(ifelse(i > 52, i - 52, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("Forecasts/Historical Average/2017-2018/EW", j, ".csv"))
  
}