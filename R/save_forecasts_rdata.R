library(tidyverse)
library(FluSight)

source("R/utils.R")
source('R/EpiDataAPI.R')

nat_1415 <- read_forecasts("Forecasts/2014-2015",
                           these_teams = c("ens-month-target-type-based-weights")) %>% 
  map_depth(1, bind_rows) %>%
  bind_rows() %>%
  mutate(season = "2014/2015")

nat_1516 <- read_forecasts("Forecasts/2015-2016",
                           these_teams = c("ens-month-target-type-based-weights")) %>% 
  map_depth(1, bind_rows) %>%
  bind_rows() %>%
  mutate(season = "2015/2016")
nat_1617 <- read_forecasts("Forecasts/2016-2017",
                           these_teams = c("ens-month-target-type-based-weights")) %>% 
  map_depth(1, bind_rows) %>%
  bind_rows() %>%
  mutate(season = "2016/2017")
nat_1718 <- read_forecasts("Forecasts/2017-2018",
                           these_teams = c("ens-month-target-type-based-weights")) %>% 
  map_depth(1, bind_rows) %>%
  bind_rows() %>%
  mutate(season = "2017/2018")
nat_1819 <- read_forecasts("Forecasts/2018-2019",
                           these_teams = c("ens-month-target-type-based-weights")) %>% 
  map_depth(1, bind_rows) %>%
  bind_rows() %>%
  mutate(season = "2018/2019")

# State forecasts
state_1415 <- read_forecasts(dir = "State Forecasts/2014-2015",
                           these_teams = c("ens-optimal-state"),
                           challenge = "state_ili") %>% 
  map_depth(1, bind_rows) %>%
  bind_rows() %>%
  mutate(season = "2014/2015")
state_1516 <- read_forecasts("State Forecasts/2015-2016",
                           these_teams = c("ens-optimal-state"),
                           challenge = "state_ili") %>% 
  map_depth(1, bind_rows) %>%
  bind_rows() %>%
  mutate(season = "2015/2016")
state_1617 <- read_forecasts("State Forecasts/2016-2017",
                           these_teams = c("ens-optimal-state"),
                           challenge = "state_ili") %>% 
  map_depth(1, bind_rows) %>%
  bind_rows() %>%
  mutate(season = "2016/2017")
state_1718 <- read_forecasts("State Forecasts/2017-2018",
                           these_teams = c("ens-optimal-state"),
                           challenge = "state_ili") %>% 
  map_depth(1, bind_rows) %>%
  bind_rows() %>%
  mutate(season = "2017/2018")
state_1819 <- read_forecasts("State Forecasts/2018-2019",
                           these_teams = c("ens-optimal-state"),
                           challenge = "state_ili") %>% 
  map_depth(1, bind_rows) %>%
  bind_rows() %>%
  mutate(season = "2018/2019")

bind_rows(nat_1415, nat_1516, nat_1617, nat_1718, nat_1819,
        state_1415, state_1516, state_1617, state_1718, state_1819) %>%
  saveRDS("Data/nat_state_forecasts.Rds")


# Create truth for each year
nat_truth_1415 <- Epidata$fluview(regions = list("nat", "hhs1", "hhs2", "hhs3",
                                               "hhs4", "hhs5", "hhs6", "hhs7",
                                               "hhs8", "hhs9", "hhs10"),
                                  epiweeks = list(Epidata$range(201440, 201528)),
                                  issue = 201528)$epidata %>%
  modify_depth(2, function(x) ifelse(is.null(x), NA, x)) %>%
  bind_rows() %>%
  mutate(week = as.integer(substr(epiweek, 5, 6)),
         location = case_when(
           region == "nat" ~ "US National",
           region == "hhs1" ~ "HHS Region 1",
           region == "hhs2" ~ "HHS Region 2",
           region == "hhs3" ~ "HHS Region 3",
           region == "hhs4" ~ "HHS Region 4",
           region == "hhs5" ~ "HHS Region 5",
           region == "hhs6" ~ "HHS Region 6",
           region == "hhs7" ~ "HHS Region 7",
           region == "hhs8" ~ "HHS Region 8",
           region == "hhs9" ~ "HHS Region 9",
           region == "hhs10" ~ "HHS Region 10"
         ),
         ili = round(ili, 1)) %>%
  select(location, week, ILI = ili) %>%
  create_truth(fluview = F, year = 2014, weekILI = .,
               start_wk = 40, end_wk = 20) %>%
  mutate(season = "2014/2015")

nat_truth_1516 <- Epidata$fluview(regions = list("nat", "hhs1", "hhs2", "hhs3",
                                               "hhs4", "hhs5", "hhs6", "hhs7",
                                               "hhs8", "hhs9", "hhs10"),
                                epiweeks = list(Epidata$range(201540, 201628)),
                                issue = 201628)$epidata %>%
  modify_depth(2, function(x) ifelse(is.null(x), NA, x)) %>%
  bind_rows() %>%
  mutate(week = as.integer(substr(epiweek, 5, 6)),
         location = case_when(
           region == "nat" ~ "US National",
           region == "hhs1" ~ "HHS Region 1",
           region == "hhs2" ~ "HHS Region 2",
           region == "hhs3" ~ "HHS Region 3",
           region == "hhs4" ~ "HHS Region 4",
           region == "hhs5" ~ "HHS Region 5",
           region == "hhs6" ~ "HHS Region 6",
           region == "hhs7" ~ "HHS Region 7",
           region == "hhs8" ~ "HHS Region 8",
           region == "hhs9" ~ "HHS Region 9",
           region == "hhs10" ~ "HHS Region 10"
         ),
         ili = round(ili, 1)) %>%
  select(location, week, ILI = ili)%>%
  create_truth(fluview = F, year = 2015, weekILI = .,
               start_wk = 40, end_wk = 20) %>%
  mutate(season = "2015/2016")

nat_truth_1617 <- Epidata$fluview(regions = list("nat", "hhs1", "hhs2", "hhs3",
                                               "hhs4", "hhs5", "hhs6", "hhs7",
                                               "hhs8", "hhs9", "hhs10"),
                                epiweeks = list(Epidata$range(201640, 201728)),
                                issue = 201728)$epidata %>%
  modify_depth(2, function(x) ifelse(is.null(x), NA, x)) %>%
  bind_rows() %>%
  mutate(week = as.integer(substr(epiweek, 5, 6)),
         location = case_when(
           region == "nat" ~ "US National",
           region == "hhs1" ~ "HHS Region 1",
           region == "hhs2" ~ "HHS Region 2",
           region == "hhs3" ~ "HHS Region 3",
           region == "hhs4" ~ "HHS Region 4",
           region == "hhs5" ~ "HHS Region 5",
           region == "hhs6" ~ "HHS Region 6",
           region == "hhs7" ~ "HHS Region 7",
           region == "hhs8" ~ "HHS Region 8",
           region == "hhs9" ~ "HHS Region 9",
           region == "hhs10" ~ "HHS Region 10"
         ),
         ili = round(ili, 1)) %>%
  select(location, week, ILI = ili) %>%
  create_truth(fluview = F, year = 2016, weekILI = .,
               start_wk = 40, end_wk = 20) %>%
  mutate(season = "2016/2017")

nat_truth_1718 <- Epidata$fluview(regions = list("nat", "hhs1", "hhs2", "hhs3",
                                               "hhs4", "hhs5", "hhs6", "hhs7",
                                               "hhs8", "hhs9", "hhs10"),
                                epiweeks = list(Epidata$range(201740, 201828)),
                                issue = 201828)$epidata %>%
  modify_depth(2, function(x) ifelse(is.null(x), NA, x)) %>%
  bind_rows() %>%
  mutate(week = as.integer(substr(epiweek, 5, 6)),
         location = case_when(
           region == "nat" ~ "US National",
           region == "hhs1" ~ "HHS Region 1",
           region == "hhs2" ~ "HHS Region 2",
           region == "hhs3" ~ "HHS Region 3",
           region == "hhs4" ~ "HHS Region 4",
           region == "hhs5" ~ "HHS Region 5",
           region == "hhs6" ~ "HHS Region 6",
           region == "hhs7" ~ "HHS Region 7",
           region == "hhs8" ~ "HHS Region 8",
           region == "hhs9" ~ "HHS Region 9",
           region == "hhs10" ~ "HHS Region 10"
         ),
         ili = round(ili, 1)) %>%
  select(location, week, ILI = ili) %>%
  create_truth(fluview = F, year = 2017, weekILI = .,
               start_wk = 40, end_wk = 20) %>%
  mutate(season = "2017/2018")

nat_truth_1819 <- Epidata$fluview(regions = list("nat", "hhs1", "hhs2", "hhs3",
                                               "hhs4", "hhs5", "hhs6", "hhs7",
                                               "hhs8", "hhs9", "hhs10"),
                                epiweeks = list(Epidata$range(201840, 201906)),
                                issue = 201906)$epidata %>%
  modify_depth(2, function(x) ifelse(is.null(x), NA, x)) %>%
  bind_rows() %>%
  mutate(week = as.integer(substr(epiweek, 5, 6)),
         location = case_when(
           region == "nat" ~ "US National",
           region == "hhs1" ~ "HHS Region 1",
           region == "hhs2" ~ "HHS Region 2",
           region == "hhs3" ~ "HHS Region 3",
           region == "hhs4" ~ "HHS Region 4",
           region == "hhs5" ~ "HHS Region 5",
           region == "hhs6" ~ "HHS Region 6",
           region == "hhs7" ~ "HHS Region 7",
           region == "hhs8" ~ "HHS Region 8",
           region == "hhs9" ~ "HHS Region 9",
           region == "hhs10" ~ "HHS Region 10"
         ),
         ili = round(ili, 1)) %>%
  select(location, week, ILI = ili) %>%
  create_truth(fluview = F, year = 2018, weekILI = .,
               start_wk = 40, end_wk = 20) %>%
  mutate(season = "2018/2019")

state_truth_1415 <- Epidata$fluview(regions = list("pa", "de", "md",
                                           "va", "wv"),
                            epiweeks = list(Epidata$range(201440, 201528)),
                            issue = 201740)$epidata %>%
  modify_depth(2, function(x) ifelse(is.null(x), NA, x)) %>%
  bind_rows() %>%
  mutate(week = as.integer(substr(epiweek, 5, 6)),
         location = case_when(
           region == "pa" ~ "Pennsylvania",
           region == "de" ~ "Delaware",
           region == "md" ~ "Maryland",
           region == "va" ~ "Virginia",
           region == "wv" ~ "West Virginia"
         ),
         ili = round(ili, 1)) %>%
  select(location, week, ILI = ili) %>%
  create_truth(fluview = F, year = 2014, weekILI = .,
               challenge = "state_ili",
               start_wk = 40, end_wk = 20) %>%
  mutate(season = "2014/2015")

state_truth_1516 <- Epidata$fluview(regions = list("pa", "de", "md",
                                           "va", "wv"),
                            epiweeks = list(Epidata$range(201540, 201628)),
                            issue = 201740)$epidata %>%
  modify_depth(2, function(x) ifelse(is.null(x), NA, x)) %>%
  bind_rows() %>%
  mutate(week = as.integer(substr(epiweek, 5, 6)),
         location = case_when(
           region == "pa" ~ "Pennsylvania",
           region == "de" ~ "Delaware",
           region == "md" ~ "Maryland",
           region == "va" ~ "Virginia",
           region == "wv" ~ "West Virginia"
         ),
         ili = round(ili, 1)) %>%
  select(location, week, ILI = ili) %>%
  create_truth(fluview = F, year = 2015, weekILI = .,
               challenge = "state_ili",
               start_wk = 40, end_wk = 20) %>%
  mutate(season = "2015/2016")

state_truth_1617 <-  Epidata$fluview(regions = list("pa", "de", "md",
                                            "va", "wv"),
                             epiweeks = list(Epidata$range(201640, 201728)),
                             issue = 201740)$epidata %>%
  modify_depth(2, function(x) ifelse(is.null(x), NA, x)) %>%
  bind_rows() %>%
  mutate(week = as.integer(substr(epiweek, 5, 6)),
         location = case_when(
           region == "pa" ~ "Pennsylvania",
           region == "de" ~ "Delaware",
           region == "md" ~ "Maryland",
           region == "va" ~ "Virginia",
           region == "wv" ~ "West Virginia"
         ),
         ili = round(ili, 1)) %>%
  select(location, week, ILI = ili) %>%
  create_truth(fluview = F, year = 2016, weekILI = .,
               challenge = "state_ili",
               start_wk = 40, end_wk = 20) %>%
  mutate(season = "2016/2017")

state_truth_1718 <-  Epidata$fluview(regions = list("pa", "de", "md",
                                            "va", "wv"),
                             epiweeks = list(Epidata$range(201740, 201828)),
                             issue = 201828)$epidata %>%
  modify_depth(2, function(x) ifelse(is.null(x), NA, x)) %>%
  bind_rows() %>%
  mutate(week = as.integer(substr(epiweek, 5, 6)),
         location = case_when(
           region == "pa" ~ "Pennsylvania",
           region == "de" ~ "Delaware",
           region == "md" ~ "Maryland",
           region == "va" ~ "Virginia",
           region == "wv" ~ "West Virginia"
         ),
         ili = round(ili, 1)) %>%
  select(location, week, ILI = ili) %>%
  create_truth(fluview = F, year = 2017, weekILI = .,
               challenge = "state_ili",
               start_wk = 40, end_wk = 20) %>%
  mutate(season = "2017/2018")

state_truth_1819 <-  Epidata$fluview(regions = list("pa", "de", "md",
                                                  "va", "wv"),
                                   epiweeks = list(Epidata$range(201840, 201906)),
                                   issue = 201906)$epidata %>%
  modify_depth(2, function(x) ifelse(is.null(x), NA, x)) %>%
  bind_rows() %>%
  mutate(week = as.integer(substr(epiweek, 5, 6)),
         location = case_when(
           region == "pa" ~ "Pennsylvania",
           region == "de" ~ "Delaware",
           region == "md" ~ "Maryland",
           region == "va" ~ "Virginia",
           region == "wv" ~ "West Virginia"
         ),
         ili = round(ili, 1)) %>%
  select(location, week, ILI = ili) %>%
  create_truth(fluview = F, year = 2018, weekILI = .,
               challenge = "state_ili",
               start_wk = 40, end_wk = 20) %>%
  mutate(season = "2018/2019")

bind_rows(nat_truth_1415, nat_truth_1516, nat_truth_1617,
          nat_truth_1718, nat_truth_1819, state_truth_1415,
          state_truth_1516, state_truth_1617, state_truth_1718,
          state_truth_1819) %>%
  saveRDS("Data/truth_from_1415.Rds")

Epidata$fluview(regions = list("nat", "hhs1", "hhs2", "hhs3",
                                          "hhs4", "hhs5", "hhs6", "hhs7",
                                          "hhs8", "hhs9", "hhs10",
                                          "pa", "de", "md", "va", "wv"),
                           epiweeks = list(Epidata$range(201840, 201906)),
                           issue = 201906)$epidata %>%
  modify_depth(2, function(x) ifelse(is.null(x), NA, x)) %>%
  bind_rows() %>%
  mutate(week = as.integer(substr(epiweek, 5, 6)),
         location = case_when(
           region == "nat" ~ "US National",
           region == "hhs1" ~ "HHS Region 1",
           region == "hhs2" ~ "HHS Region 2",
           region == "hhs3" ~ "HHS Region 3",
           region == "hhs4" ~ "HHS Region 4",
           region == "hhs5" ~ "HHS Region 5",
           region == "hhs6" ~ "HHS Region 6",
           region == "hhs7" ~ "HHS Region 7",
           region == "hhs8" ~ "HHS Region 8",
           region == "hhs9" ~ "HHS Region 9",
           region == "hhs10" ~ "HHS Region 10",
           region == "pa" ~ "Pennsylvania",
           region == "de" ~ "Delaware",
           region == "md" ~ "Maryland",
           region == "va" ~ "Virginia",
           region == "wv" ~ "West Virginia"
         ),
         ili = round(ili, 1)) %>%
  select(location, week, ILI = ili) %>%
  saveRDS("Data/obs_ili_1819.Rds")
