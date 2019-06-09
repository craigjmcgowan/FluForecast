# Generate subtype-weighted historical average forecast
library(tidyverse)
library(FluSight)
library(zoo)

# Load data -------
load("Data/ili_state.Rdata")
load("Data/virologic_state.Rdata")
load("Data/virologic.Rdata")

# Load functions 
source("R/create_subtype_forecast.R")
source("R/utils.R")

# Calculate cumulative virologic status across seasons -----
vir_ssn_per <- state_virologic %>%
  group_by(season, location) %>%
  slice(tail(row_number(), 1)) %>%
  select(season, location, h1per = h1per_of_a, h3per = h3per_of_a) %>%
  ungroup() %>%
  filter(!location %in% c("New Jersey", "Rhode Island")) %>%
  bind_rows(
    # For NJ and RI - bring in HHS region percentages to use instead
    virologic_combined %>%
      filter(location %in% c("HHS Region 1", "HHS Region 2"),
             season %in% c("2010/2011", "2011/2012", "2012/2013", "2013/2014",
                           "2014/2015", "2015/2016", "2016/2017", "2017/2018",
                           "2018/2019")) %>%
      group_by(season, location) %>%
      slice(tail(row_number(), 1)) %>%
      select(season, location, h1per = h1per_of_a, h3per = h3per_of_a) %>%
      ungroup() %>%
      mutate(location = case_when(location == "HHS Region 1" ~ "Rhode Island",
                                  location == "HHS Region 2" ~ "New Jersey"))
  )


# Create forecasts for 2014/2015 ------
train_ili_1415 <- state_current %>%
  filter(year <= 2014, season != "2014/2015")

virologic_1415 <- state_virologic %>%
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
  select(year,location, week, h1sum, h1per, h3sum, h3per, order_week) %>%
  ungroup()

# Create directory to store forecasts
dir.create("State Forecasts/2014-2015/Subtype Historical Average",
           showWarnings = FALSE)

# Create target densities and functions
subtype_densities_1415 <- create_subtype_densities(
  train_ili_1415, vir_ssn_per, challenge = "state_ili"
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
    challenge = "state_ili"
  )
  
  j <- str_pad(ifelse(i > 53, i - 53, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("State Forecasts/2014-2015/Subtype Historical Average/EW", j, ".csv"))
  
}

# Create forecasts for 2015/2016 ------
train_ili_1516 <- state_current %>%
  filter(year <= 2015, season != "2015/2016")

virologic_1516 <- state_virologic %>%
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
  select(year,location, week, h1sum, h1per, h3sum, h3per, order_week) %>%
  ungroup()

# Create directory to store forecasts
dir.create("State Forecasts/2015-2016/Subtype Historical Average",
           showWarnings = FALSE)

# Create target densities and functions
subtype_densities_1516 <- create_subtype_densities(
  train_ili_1516, vir_ssn_per, challenge = "state_ili"
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
    challenge = "state_ili"
  )
  
  j <- str_pad(ifelse(i > 52, i - 52, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("State Forecasts/2015-2016/Subtype Historical Average/EW", j, ".csv"))
  
}

# Create forecasts for 2016/2017 ------
train_ili_1617 <- state_current %>%
  filter(year <= 2016, season != "2016/2017")

virologic_1617 <- state_virologic %>%
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
  select(year,location, week, h1sum, h1per, h3sum, h3per, order_week) %>%
  ungroup()

# Create directory to store forecasts
dir.create("State Forecasts/2016-2017/Subtype Historical Average",
           showWarnings = FALSE)

# Create target densities and functions
subtype_densities_1617 <- create_subtype_densities(
  ili_df = train_ili_1617, vir_ssn_per = vir_ssn_per, 
  challenge = "state_ili"
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
    challenge = "state_ili"
  )
  
  j <- str_pad(ifelse(i > 52, i - 52, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("State Forecasts/2016-2017/Subtype Historical Average/EW", j, ".csv"))
  
}

# Create forecasts for 2017/2018 ------
train_ili_1718 <- state_current %>%
  filter(year <= 2017, season != "2017/2018",
         location != "Louisiana") # Remove LA b/c only one season of training data

virologic_1718 <- state_virologic %>%
  filter(season == "2017/2018", location != "Louisiana") %>%
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
  select(year,location, week, h1sum, h1per, h3sum, h3per, order_week) %>%
  ungroup()

# Create directory to store forecasts
dir.create("State Forecasts/2017-2018/Subtype Historical Average",
           showWarnings = FALSE)

# Create target densities and functions
subtype_densities_1718 <- create_subtype_densities(
  train_ili_1718, vir_ssn_per, challenge = "state_ili"
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
    challenge = "state_ili"
  )
  
  j <- str_pad(ifelse(i > 52, i - 52, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("State Forecasts/2017-2018/Subtype Historical Average/EW", j, ".csv"))
  
}

# Create forecasts for 2018/2019 ------
train_ili_1819 <- state_current %>%
  filter(year <= 2018, season != "2018/2019")

virologic_1819 <- state_virologic %>%
  filter(season == "2018/2019") %>%
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
  select(year,location, week, h1sum, h1per, h3sum, h3per, order_week) %>%
  ungroup()

# Create directory to store forecasts
dir.create("State Forecasts/2018-2019/Subtype Historical Average",
           showWarnings = FALSE)

# Create target densities and functions
subtype_densities_1819 <- create_subtype_densities(
  train_ili_1819, vir_ssn_per, challenge = "state_ili"
)
subtype_functions_1819 <- modify_depth(
  subtype_densities_1819, 3,
  function(dens) approxfun(dens$x, dens$y, rule = 2)
)
MMWRweek(Sys.Date())
for(i in 65:69) {
  temp <- create_subtype_forecast(
    functions = subtype_functions_1819,
    virologic = virologic_1819,
    pub_week = i,
    season = "2018/2019",
    challenge = "state_ili"
  )
  
  j <- str_pad(ifelse(i > 52, i - 52, i), 2, pad = "0")
  
  write_csv(temp,
            path = paste0("State Forecasts/2018-2019/Subtype Historical Average/EW", j, ".csv"))
  
}
