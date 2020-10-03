# Generate scores
library(tidyverse)
library(FluSight)

source("R/utils.R")

# Helper function to create truth from Delphi-formatted file
read_ili <- function(pub_week, pub_season) {
  
  readRDS("Data/ili_init_pub_list.rds")[[pub_week]] %>%
    filter(season == pub_season, !location %in% state.name, !location %in% c('District of Columbia', 'Puerto Rico')) %>%
    mutate(order_week = week_inorder(week, pub_season)) %>%
    arrange(location, order_week) %>%
    select(location, week, ILI) %>%
    mutate(ILI = round(ILI, 1))
  
}

truth_epidata <- function(pub_week, pub_season) {
  
  create_truth(fluview = FALSE, year = as.numeric(substr(pub_season, 1, 4)), 
              weekILI = read_ili(pub_week, pub_season), challenge = "ilinet", 
              start_wk = 40, end_wk = 20)
  
}

restrict_forecasts <- function(forecasts, truth_week) {
  
  weeks <- c("EW40", "EW41", "EW42", "EW43", "EW44", "EW45", "EW46", "EW47", 
             "EW48", "EW49", "EW50", "EW51", "EW52", "EW01", "EW02", "EW03", 
             "EW04", "EW05", "EW06", "EW07", "EW08", "EW09", "EW10", "EW11", 
             "EW12", "EW13", "EW14", "EW15", "EW16", "EW17", "EW18", "EW19", 
             "EW20")
  
  teams <- names(forecasts)
  if(truth_week %in% c('EW21', 'EW22')) {
    valid_weeks <- weeks
  } else {
    valid_weeks <- weeks[1:which(weeks == truth_week) - 1]
  }
  
  out_list <- list()
  
  for(this_team in teams) {
    for(this_week in valid_weeks) {
      out_list[[this_team]][[this_week]] <- forecasts[[this_team]][[this_week]]
    }
  }
  
  return(out_list)
}

# Read in retrospective forecasts -------
forecasts_1819 <- read_forecasts("Forecasts/Live/2018-2019")
forecasts_1920 <- read_forecasts("Forecasts/Live/2019-2020")

# Create weekly truths
ili_1819 <- map(weeks_1819$epidata_week, read_ili, '2018/2019')
names(ili_1819) <- weeks_1819$flusight_week
ili_1920 <- map(weeks_1920$epidata_week, read_ili, '2019/2020')
names(ili_1920) <- weeks_1920$flusight_week

weeks_1819 <- tibble(
  epidata_week = c(as.character(seq(201841, 201852)), as.character(seq(201901, 201922))),
  flusight_week = paste0('EW', substr(epidata_week, 5, 6))
)
weeks_1920 <- tibble(
  epidata_week = c(as.character(seq(201941, 201952)), as.character(seq(202001, 202022))),
  flusight_week = paste0('EW', substr(epidata_week, 5, 6))
)

truth_1819 <- map(weeks_1819$epidata_week, truth_epidata, '2018/2019')
names(truth_1819) <- weeks_1819$flusight_week
truth_1920 <- map(weeks_1920$epidata_week, truth_epidata, '2019/2020')
names(truth_1920) <- weeks_1920$flusight_week

# Calculate model scores as known each week
weekly_scores_1819 <- list()
weekly_eval_scores_1819 <- list()
weekly_scores_1920 <- list()
weekly_eval_scores_1920 <- list()

for (i in seq_along(truth_1819)) {
  
  week_name <- names(truth_1819[i])
  week_data <- truth_1819[[i]]
  ili_data <- ili_1819[[i]]
  
  eval_period <- create_eval_period(ili_data, week_data, '2018/2019')
  
  observed_truth <- filter(week_data, !is.na(bin_start_incl), bin_start_incl != 'none')

  forecasts_to_score <- restrict_forecasts(forecasts_1819, week_name)
  
  weekly_scores_1819[[week_name]] <- calc_scores(forecasts_to_score, observed_truth, '2018/2019', exclude = T)
  weekly_eval_scores_1819[[week_name]] <- create_eval_scores(weekly_scores_1819[[week_name]], eval_period, '2018/2019')
  
}

for (i in seq_along(truth_1920)) {
  
  week_name <- names(truth_1920[i])
  week_data <- truth_1920[[i]]
  ili_data <- ili_1920[[i]]
  
  eval_period <- create_eval_period(ili_data, week_data, '2019/2020')
  
  observed_truth <- filter(week_data, !is.na(bin_start_incl), bin_start_incl != 'none')
  
  forecasts_to_score <- restrict_forecasts(forecasts_1920, week_name)
  
  weekly_scores_1920[[week_name]] <- calc_scores(forecasts_to_score, observed_truth, '2019/2020', exclude = T)
  weekly_eval_scores_1920[[week_name]] <- create_eval_scores(weekly_scores_1920[[week_name]], eval_period, '2019/2020')
  
}

saveRDS(weekly_scores_1819, 'Data/weekly_scores_1819.RDS')
saveRDS(weekly_eval_scores_1819, 'Data/weekly_eval_scores_1819.RDS')
saveRDS(weekly_scores_1920, 'Data/weekly_scores_1920.RDS')
saveRDS(weekly_eval_scores_1920, 'Data/weekly_eval_scores_1920.RDS')
