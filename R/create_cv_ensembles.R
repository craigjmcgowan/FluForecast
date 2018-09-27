# Generate ensemble forecasts
library(tidyverse)
library(FluSight)

source("R/utils.R")

# List forecast files to be weighted
model_files <- list.files(path = "Forecasts", recursive = TRUE, full.names = T)
model_files <- model_files[!grepl('ens', model_files)]

# List weighting schemes to be calculated
weight_files <- list.files("Weights")
weight_types <- sub(
  pattern = "-weights.csv", 
  replacement = "",
  weight_files)

for(j in 1:length(weight_files)) {
  stacking_weights <- read.csv(paste0("weights/", weight_files[j]), 
                               stringsAsFactors=FALSE)
  stacked_name <- sub(pattern = ".csv", replacement = "", weight_files[j])
  
  # Remove future season from weights
  seasons <- unique(stacking_weights$season)
  if("2018/2019" %in% seasons)
    seasons <- seasons[-which(seasons=="2018/2019")]
  
  # If weights by week included, reset to MMWR week
  if (any(grepl("Model.Week", names(stacking_weights)))) {
    stacking_weights <- mutate(
      stacking_weights,
      `Model.Week` = week_reset(`Model.Week`, season),
      ew = paste0("EW", str_pad(`Model.Week`, 2, "left", 0))
    ) %>%
      select(-`Model.Week`) %>%
      # Ad hoc correction to make weights for EW20 in 2014/2015
      bind_rows(., filter(., season == "2014/2015",
                          ew == "EW19") %>%
                  mutate(ew = "EW20"))
  }
  
  ## loop through each season and each season-week to make stacked forecasts
  for(i in 1:length(seasons)){
    
    loso_season =  seasons[i]
    wt_subset <- filter(stacking_weights, season==loso_season) %>%
      dplyr::select(-season)
    
    dir.create(file.path("Forecasts", sub("/", "-", loso_season),
                         paste0("ens-", stacked_name)), 
               showWarnings = FALSE)
    
    ## identify the "EWXX-YYYY" combos for files given season
    first_year <- substr(loso_season, 0, 4)
    first_year_season_weeks <- if(first_year==2014) {40:53} else {40:52}
    week_names <- c(paste0("EW", first_year_season_weeks),
                    paste0("EW", str_pad(1:20, 2, "left", pad = 0)))

    for (k in 1:length(week_names)) {
  
      this_week <- week_names[k]
      message(paste(stacked_name, "::", this_week, "::", Sys.time()))
      
      if (any(grepl("ew", names(wt_subset)))) {
        wt_sub_subset <- filter(wt_subset, ew == this_week) %>%
          select(-ew)
      } else {
        wt_sub_subset <- wt_subset
      }
      
      ## stack models, save ensemble file
      files_to_stack <- model_files[(grepl(this_week, model_files) & 
                                      grepl(sub("/", "-", loso_season), model_files))]

      file_df <- data.frame(
        file = files_to_stack, 
        model_id = word(files_to_stack, start = 3, end = 3, sep = "/"),
        stringsAsFactors = FALSE)
      
      stacked_entry <- stack_forecasts(file_df, wt_sub_subset)
      stacked_file_name <- paste0(
        "Forecasts/", sub("/", "-", loso_season), "/ens-",
        stacked_name, "/", this_week, ".csv"
      )
      write.csv(stacked_entry, file=stacked_file_name, 
                row.names = FALSE, quote = FALSE)
    }
  }
}

