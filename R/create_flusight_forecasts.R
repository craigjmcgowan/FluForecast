# Create past forecasts for FluSight Network
library(tidyverse)
library(forecast)
library(lubridate)
library(FluSight)
library(MMWRweek)

# Load functions
source("R/utils.R")

##### Load data ######
load("Data/ili.Rdata")
load("Data/virologic.Rdata")
load("Data/Gtrends.Rdata")
load("Data/state_gtrend.Rdata")

##### Create truth for all seasons and combine datasets #####
load("Data/truth_and_data.Rdata")
load("Data/CV_Transform_terms.Rdata")
load("data/CV_Fourier_terms.Rdata")
load("Data/CV_ARIMA_terms.Rdata")
load("Data/CV_covar_terms.Rdata")

##### Subtype-weighted historical average model (Kudu) #####
for(this_season in c("2010-2011", "2011-2012", "2012-2013", "2013-2014",
                     "2014-2015", "2015-2016", "2016-2017", "2017-2018")) {
  subtype_files <- list.files(paste0("Forecasts/", this_season,
                                     "/Subtype Historical Average"))

  for(i in seq_along(subtype_files)) {
    
    print_week <- sub(pattern = ".csv", replacement = "", subtype_files[i])
    
    print_year <- ifelse(as.numeric(substr(print_week, 3, 4)) < 40,
                         substr(this_season, 6, 9),
                         substr(this_season, 1, 4))
    
    read_csv(paste0("Forecasts/", this_season, "/Subtype Historical Average/", 
                    subtype_files[i])) %>%
      select(location, target, type, unit, bin_start_incl, bin_end_notincl, value) %>%
      write_csv(paste0("../cdc-flusight-ensemble/model-forecasts/component-models/Protea_Kudu/",
                       print_week, "-", print_year, "-Protea_Kudu.csv"))
  }
  
}

##### Dynamic Harmonic Model (Springbok) #####
load("Data/Final_CV_fits.Rdata")

final_forecast_data_setup <- crossing(season = c("2010/2011", "2011/2012", "2012/2013",
                                                 "2013/2014", "2014/2015", "2015/2016",
                                                 "2016/2017", "2017/2018"),
                                      week = c(40:73),
                                      location = unique(flu_data_merge$location)) %>%
  filter(week < 73 | season == "2014/2015") %>%
  mutate(
    epiweek = case_when(
      season == "2014/2015" & week > 53 ~ 
        as.numeric(paste0(substr(season, 6, 9), 
                          str_pad(week - 53, 2, "left", "0"))),
      week > 52 ~ 
        as.numeric(paste0(substr(season, 6, 9), 
                          str_pad(week - 52, 2, "left", "0"))),
      TRUE ~ as.numeric(paste0(substr(season, 1, 4), str_pad(week, 2, "left", "0")))
    ),
    prev_season = case_when(
      season == "2010/2011" ~ "2007/2008",
      TRUE ~ paste0(as.numeric(substr(season, 1, 4)) - 1, "/",
                    substr(season, 1, 4))
    )
  ) %>%
  inner_join(final_fits, 
             by = c("location", "season")) 

springbok_path <- "../cdc-flusight-ensemble/model-forecasts/component-models/Protea_Springbok" 

# Create and save forecast files
for(this_season in c("2016/2017", "2017/2018")) { #unique(final_forecast_data_setup$season)) {
  
  for(this_week in unique(final_forecast_data_setup$week[final_forecast_data_setup$season == this_season])) {

    temp <- final_forecast_data_setup %>%
      filter(season == this_season, week == this_week) %>%
      mutate(
        pred_data = pmap(list(season, week, location, epiweek), 
                         ~ filter(flu_data_merge, year <= as.numeric(substr(..1, 6, 9)),
                                  season != paste0(substr(..1, 6, 9), "/",
                                                   as.numeric(substr(..1, 6, 9)) + 1),
                                  season != ..1 | order_week %in% 40:..2,
                                  location == ..3) %>%
                           left_join(select(ili_init_pub_list[[paste(..4)]],
                                            ILI, epiweek, location),
                                     by = c("epiweek", "location")) %>%
                           mutate(ILI = ifelse(is.na(ILI.y), ILI.x, ILI.y)) %>%
                           select(-ILI.x, -ILI.y) %>%
                           mutate(ILI = ts(ILI, frequency = 52, start = c(2006, 40)))),
        # Set up Fourier data for forecasting
        xreg = case_when(
          model == "ARIMA only" ~ 
            map2(pred_data, K,
                 ~ as.data.frame(fourier(.x$ILI, K = .y))),
          model == "National Gtrends" ~
            pmap(list(pred_data, K), 
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2)),
                         ..1$hits)),
          model == "Regional Gtrends" ~
            pmap(list(pred_data, K), 
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2)),
                         ..1$region_hits)),
          model == "FluVirus" ~ 
            pmap(list(pred_data, K),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2)),
                         ..1$cum_h1per, ..1$cum_h3per)),
          model == "Backfill" ~ 
            pmap(list(pred_data, K),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2)),
                         ..1$backfill)),
          model == "National Gtrends & Backfill" ~
            pmap(list(pred_data, K), 
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2)),
                         ..1$hits, ..1$backfill)),
          model == "Regional Gtrends & Backfill" ~
            pmap(list(pred_data, K), 
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2)),
                         ..1$region_hits, ..1$backfill)),
          model == "FluVirus & Backfill" ~ 
            pmap(list(pred_data, K),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2)),
                         ..1$cum_h1per, ..1$cum_h3per,
                         ..1$backfill)),
          model == "National Gtrends & FluVirus" ~ 
            pmap(list(pred_data, K),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2)),
                         ..1$hits, ..1$cum_h1per, 
                         ..1$cum_h3per)),
          model == "Regional Gtrends & FluVirus" ~ 
            pmap(list(pred_data, K),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2)),
                         ..1$region_hits, ..1$cum_h1per, 
                         ..1$cum_h3per)),
          model == "National Gtrends, FluVirus, & Backfill" ~ 
            pmap(list(pred_data, K),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2)),
                         ..1$hits, ..1$cum_h1per, 
                         ..1$cum_h3per, ..1$backfill)),
          model == "Regional Gtrends, FluVirus, & Backfill" ~ 
            pmap(list(pred_data, K),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2)),
                         ..1$region_hits, ..1$cum_h1per, 
                         ..1$cum_h3per, ..1$backfill))
        ),
        max_week = ifelse(season == "2014/2015", 53, 52),
        # Fit models
        pred_fit = pmap(list(pred_data, fit, xreg),
                        ~ Arima(..1$ILI, xreg = ..3, model = ..2)),
        # Create data frame of xreg terms for forecasting
        gtrend_forecast = pmap(
          list(pred_data, location, season, prev_season, week, max_week),
          ~ tibble(hits = c(flu_data_merge %>%
                              filter(location == ..2, season == ..3, 
                                     order_week == ..5 + 1) %>%
                              pull(hits),
                            flu_data_merge %>%
                              filter(location == ..2, season == ..4, 
                                     order_week > ..5 + 1, 
                                     order_week < ..6 + 26) %>%
                              pull(hits) *
                              mean(flu_data_merge %>%
                                     filter(location == ..2, season == ..3,
                                            order_week <= ..5 + 1) %>%
                                     pull(hits) /
                                     flu_data_merge %>%
                                     filter(location == ..2, season == ..4,
                                            order_week <= ..5 + 1) %>%
                                     pull(hits))))
        ),
        reg_gtrend_forecast = pmap(
          list(pred_data, location, season, prev_season, week, max_week),
          ~ tibble(region_hits = c(flu_data_merge %>%
                                     filter(location == ..2, season == ..3, 
                                            order_week == ..5 + 1) %>%
                                     pull(region_hits),
                                   flu_data_merge %>%
                                     filter(location == ..2, season == ..4, 
                                            order_week > ..5 + 1, 
                                            order_week < ..6 + 26) %>%
                                     pull(region_hits) *
                                     mean(flu_data_merge %>%
                                            filter(location == ..2, season == ..3,
                                                   order_week <= ..5 + 1) %>%
                                            pull(region_hits) /
                                            flu_data_merge %>%
                                            filter(location == ..2, season == ..4,
                                                   order_week <= ..5 + 1) %>%
                                            pull(region_hits))))
        ),
        h1_per_forecast = pmap(
          list(pred_data, season, max_week),
          ~ rep(last(..1$cum_h1per), ..3 - 14 - 
                  nrow(..1[..1$season == ..2, ]))
        ),
        h3_per_forecast = pmap(
          list(pred_data, season, max_week),
          ~ rep(last(..1$cum_h3per), ..3 - 14 - 
                  nrow(..1[..1$season == ..2, ]))
        ),
        backfill = pmap(
          list(data, week, max_week, season),
          function(data, week, max_week, season) {
            out <- numeric()
            for(i in week:(max_week + 24)) {
              
              if(season %in% c("2010/2011", "2011/2012", "2012/2013")) {
                temp_data <- data[data$order_week == i, ]
              } else{
                temp_data <- data[data$order_week == i &
                               !data$season1 %in% c("2004/2005", "2005/2006", 
                                                   "2006/2007", "2007/2008"), ]
              }
              
              out <- c(out, predict(lm(backfill ~ year, 
                                       data = temp_data),
                                    data.frame(year = as.numeric(substr(season, 1, 4)))))
            }
            out
          } 
        ),
        forecast_xreg = case_when(
          model == "ARIMA only" ~ 
            pmap(list(pred_data, K, max_week, season),
                 ~ as.data.frame(fourier(..1$ILI, K = ..2,
                                         h = ..3 - 14 - 
                                           nrow(..1[..1$season == ..4, ])))),
          model == "Backfill" ~
            pmap(list(pred_data, K, max_week, season, backfill),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                               h = ..3 - 14 - 
                                                 nrow(..1[..1$season == ..4, ]))),
                         data.frame(backfill = ..5)) %>%
                   rename(`..1$backfill` = backfill)),
          model == "National Gtrends" ~
            pmap(list(pred_data, K, max_week, season, gtrend_forecast),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                               h = ..3 - 14 - 
                                                 nrow(..1[..1$season == ..4, ]))),
                         data.frame(hits = ..5)) %>%
                   rename(`..1$hits` = hits)),
          model == "Regional Gtrends" ~
            pmap(list(pred_data, K, max_week, season, reg_gtrend_forecast),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                               h = ..3 - 14 - 
                                                 nrow(..1[..1$season == ..4, ]))),
                         data.frame(region_hits = ..5)) %>%
                   rename(`..1$region_hits` = region_hits)),
          model == "FluVirus" ~
            pmap(list(pred_data, K, max_week, season, 
                      h1_per_forecast, h3_per_forecast),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                               h = ..3 - 14 - 
                                                 nrow(..1[..1$season == ..4, ]))),
                         data.frame(h1_per_samples = ..5,
                                    h3_per_samples = ..6)) %>%
                   rename(`..1$cum_h1per` = h1_per_samples,
                          `..1$cum_h3per` = h3_per_samples)),
          model == "National Gtrends & Backfill" ~
            pmap(list(pred_data, K, max_week, season, gtrend_forecast, backfill),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                               h = ..3 - 14 - 
                                                 nrow(..1[..1$season == ..4, ]))),
                         data.frame(hits = ..5,
                                    backfill = ..6)) %>%
                   rename(`..1$hits` = hits,
                          `..1$backfill` = backfill)),
          
          model == "Regional Gtrends & Backfill" ~
            pmap(list(pred_data, K, max_week, season, reg_gtrend_forecast, backfill),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                               h = ..3 - 14 - 
                                                 nrow(..1[..1$season == ..4, ]))),
                         data.frame(region_hits = ..5,
                                    backfill = ..6)) %>%
                   rename(`..1$region_hits` = region_hits,
                          `..1$backfill` = backfill)),
          model == "FluVirus & Backfill" ~
            pmap(list(pred_data, K, max_week, season, 
                      h1_per_forecast, h3_per_forecast, backfill),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                               h = ..3 - 14 - 
                                                 nrow(..1[..1$season == ..4, ]))),
                         data.frame(h1_per_samples = ..5,
                                    h3_per_samples = ..6,
                                    backfill = ..7)) %>%
                   rename(`..1$cum_h1per` = h1_per_samples,
                          `..1$cum_h3per` = h3_per_samples,
                          `..1$backfill` = backfill)),
          model == "National Gtrends & FluVirus" ~
            pmap(list(pred_data, K, max_week, season, gtrend_forecast,
                      h1_per_forecast, h3_per_forecast),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                               h = ..3 - 14 - 
                                                 nrow(..1[..1$season == ..4, ]))),
                         data.frame(hits = ..5, 
                                    h1_per_samples = ..6,
                                    h3_per_samples = ..7)) %>%
                   rename(`..1$hits` = hits,
                          `..1$cum_h1per` = h1_per_samples,
                          `..1$cum_h3per` = h3_per_samples)),
          model == "Regional Gtrends & FluVirus" ~
            pmap(list(pred_data, K, max_week, season, reg_gtrend_forecast,
                      h1_per_forecast, h3_per_forecast),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                               h = ..3 - 14 - 
                                                 nrow(..1[..1$season == ..4, ]))),
                         data.frame(region_hits = ..5, 
                                    h1_per_samples = ..6,
                                    h3_per_samples = ..7)) %>%
                   rename(`..1$region_hits` = region_hits,
                          `..1$cum_h1per` = h1_per_samples,
                          `..1$cum_h3per` = h3_per_samples)),
          model == "National Gtrends, FluVirus, & Backfill" ~
            pmap(list(pred_data, K, max_week, season, gtrend_forecast,
                      h1_per_forecast, h3_per_forecast, backfill),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                               h = ..3 - 14 - 
                                                 nrow(..1[..1$season == ..4, ]))),
                         data.frame(hits = ..5, 
                                    h1_per_samples = ..6,
                                    h3_per_samples = ..7,
                                    backfill = ..8)) %>%
                   rename(`..1$hits` = hits,
                          `..1$cum_h1per` = h1_per_samples,
                          `..1$cum_h3per` = h3_per_samples,
                          `..1$backfill` = backfill)),
          model == "Regional Gtrends, FluVirus, & Backfill" ~
            pmap(list(pred_data, K, max_week, season, reg_gtrend_forecast,
                      h1_per_forecast, h3_per_forecast, backfill),
                 ~ cbind(as.data.frame(fourier(..1$ILI, K = ..2,
                                               h = ..3 - 14 - 
                                                 nrow(..1[..1$season == ..4, ]))),
                         data.frame(region_hits = ..5, 
                                    h1_per_samples = ..6,
                                    h3_per_samples = ..7, 
                                    backfill = ..8)) %>%
                   rename(`..1$region_hits` = region_hits,
                          `..1$cum_h1per` = h1_per_samples,
                          `..1$cum_h3per` = h3_per_samples,
                          `..1$backfill` = backfill))
        ),
        pred_results = pmap(
          list(pred_fit, pred_data, forecast_xreg, location, season, max_week),
          ~ fit_to_forecast(object = ..1,
                            xreg = ..3,
                            pred_data = ..2,
                            location = ..4,
                            season = ..5,
                            max_week = ..6,
                            npaths = 1000)
        )
      ) %>%
      select(location, week, max_week, pred_results)
    
    EW <- case_when(temp$week[1] > temp$max_week[1] ~ 
                      str_pad(temp$week[1] - temp$max_week[1], 2, side = "left", 0),
                    TRUE ~ str_pad(temp$week[1], 2, "left", 0))
    
    this_year <- case_when(
      as.numeric(EW) >= 40 ~ substr(this_season, 1, 4),
      as.numeric(EW) < 40 ~ substr(this_season, 6, 9)
    )
    
    
    temp %>%
      unnest() %>%
      select(-week, -max_week, -forecast_week) %>%
      bind_rows(generate_point_forecasts(.)) %>%
      select(location, target, type, unit, bin_start_incl, bin_end_notincl, value) %>%
      write_csv(path = paste0(
        springbok_path, "/EW", EW, "-", this_year, "-Protea_Springbok.csv"))
    
  }
}





##### Weighted Ensemble Model (Cheetah) #####

# List forecast files to be weighted
model_files <- list.files(path = "Forecasts", recursive = TRUE, full.names = T)
model_files <- model_files[!grepl('ens', model_files)]

# List weighting schemes to be calculated
weight_files <- "Weights/month-target-type-based-weights.csv"
weight_types <- sub(
  pattern = "-weights.csv", 
  replacement = "",
  sub(pattern = "Weights/",
      replacement = "",
      weight_files)
)

# Read in weights
stacking_weights <- read.csv(weight_files, stringsAsFactors=FALSE)
stacked_name <- sub(pattern = ".csv", replacement = "", weight_files)

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

sub_folder <- "../cdc-flusight-ensemble/model-forecasts/component-models/Protea_Cheetah"

dir.create(sub_folder, showWarnings = FALSE)

## loop through each season and each season-week to make stacked forecasts
for(i in 1:length(seasons)){

  loso_season =  seasons[i]
  wt_subset <- filter(stacking_weights, season==loso_season) %>%
    dplyr::select(-season)
  
  ## identify the "EWXX-YYYY" combos for files given season
  first_year <- substr(loso_season, 0, 4)
  second_year <- substr(loso_season, 6, 9)
  first_year_season_weeks <- if(first_year==2014) {40:53} else {40:52}
  week_names <- c(paste0("EW", first_year_season_weeks),
                  paste0("EW", str_pad(1:20, 2, "left", pad = 0)))
  
  for (k in 1:length(week_names)) {
   
    this_week <- week_names[k]
    this_year <- case_when(
      as.numeric(substr(week_names[k], 3, 4)) >= 40 ~ first_year,
      as.numeric(substr(week_names[k], 3, 4)) < 40 ~ second_year
    )
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
    
    stacked_entry <- stack_forecasts(file_df, wt_sub_subset) %>%
      select(location, target, type, unit, bin_start_incl, bin_end_notincl, value)
    
    stacked_file_name <- paste0(
      sub_folder, "/", this_week, "-", this_year, "-Protea_Cheetah.csv"
    )
    
    write.csv(stacked_entry, file=stacked_file_name, 
              row.names = FALSE, quote = FALSE)
  }
}

