# Create weekly prospective forecasts for 2019-2020 season

library(tidyverse)
library(FluSight)
library(MMWRweek)
library(cdcfluview)
library(forecast)
library(gtrendsR)
library(lubridate)
library(zoo)

source("R/utils.R")
source("R/create_subtype_forecast.R")

##### Set week that forecasts are being based on #####
EW <- 47
epiweek <- 201947
EW_paste <- str_pad(EW, 2, pad = "0")
order_week <- ifelse(EW < 40, EW + 52, EW)

##### Update data #####
source("R/read_data.R")
source("R/save_who_nrevss.R")

# Remove state info
ili_current <- filter(ili_current, !location %in% state.name)
virologic_combined <- filter(virologic_combined, !location %in% state.name)
gtrend <- gtrend %>%
  # Create HHS region proxies for Gtrend data
  mutate(location = case_when(
    location == "Massachusetts" ~ "HHS Region 1",
    location == "New York" ~ "HHS Region 2",
    location == "Pennsylvania" ~ "HHS Region 3",
    location == "Florida" ~ "HHS Region 4",
    location == "Illinois" ~ "HHS Region 5",
    location == "Texas" ~ "HHS Region 6",
    location == "Missouri" ~ "HHS Region 7",
    location == "Colorado" ~ "HHS Region 8",
    location == "California" ~ "HHS Region 9",
    location == "Washington" ~ "HHS Region 10",
    TRUE ~ location
  )) %>%
  filter(!location %in% state.name, !season %in% c("2008/2009", "2009/2010"),
         !(year == 2015 & week == 33))

# Load backfill densities
ili_backfill_densities <- readRDS('Data/ili_backfill_densities.Rds') %>%
  filter(!location %in% state.name, season == '2019/2020')

flu_data_merge <- select(ili_current, epiweek, ILI, year, week, season, location) %>%
  # Add virologic data
  inner_join(select(virologic_combined, location, season, year, week, cum_h1per_6wk,
                    cum_h3per_6wk, cum_bper_6wk),
             by = c("location", "season", "year", "week")) %>%
  # Add Google Trends data by location
  right_join(select(gtrend, season, week, year, location, region_hits = hits),
             by = c("location", "season", "year", "week")) %>%
  right_join(filter(gtrend, location == "US National") %>%
               select(season, week, year, hits),
             by = c("season", "year", "week")) %>%
  # Remove 2008/2009 and 2009/2010 seasons due to pandemic activity
  filter(!season %in% c("2008/2009", "2009/2010")) %>%
  # Remove week 33 in 2015 so all seasons have 52 weeks - minimal activity
  filter(!(year == 2015 & week == 33)) %>%
  mutate(order_week = week_inorder(week, season)) 
  
##### Steenbok ######

# Move current week's forecast to Shiny app folder
steenbok_pred <- read_entry(paste0('Forecasts/Live/2019-2020/Historical Average/EW',
                                   EW_paste, '.csv')) %>%
  select(-forecast_week)

write_csv(steenbok_pred, 
          paste0('../ForecastShiny/RawData/2019-2020/NatReg/Historical Average/EW',
                 EW_paste, '.csv'))

##### Kudu #####
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

# Create data.frame of probability of no onset based on prior years 
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


# Create forecasts 
kudu_ili <- ili_current %>%
  filter(year <= 2019, season != "2019/2020")

kudu_virologic <- virologic_combined %>%
  filter(season == "2019/2020", !location %in% state.name) %>%
  # Create variable to keep weeks in order
  mutate(order_week = week_inorder(week, season)) %>%
  select(year, location, week, order_week, 
         h1per_wk = cum_h1per_6wk, h3per_wk = cum_h3per_6wk,
         bper_wk = cum_bper_6wk, h1per_ssn = cum_h1per,
         h3per_ssn = cum_h3per, bper_ssn = cum_bper)

# Create directory to store forecasts
dir.create("Forecasts/Live/2019-2020/Subtype Historical Average/",
           showWarnings = FALSE)

kudu_flusight_path <- 
  "../cdc-flusight-ensemble/model-forecasts/real-time-component-models/Protea_Kudu" 

kudu_shiny_path <- 
  '../ForecastShiny/RawData/2019-2020/NatReg/Subtype Historical Average'

# Create target densities and functions
subtype_densities_1920 <- create_subtype_densities(
  ili_df = kudu_ili, vir_ssn_per = vir_ssn_per, vir_wk_per = vir_wk_per
)

subtype_functions_1920 <- modify_depth(
  subtype_densities_1920, 3,
  function(dens) approxfun(dens$x, dens$y, rule = 2)
)

kudu_pred <- create_subtype_forecast(
    functions = subtype_functions_1920,
    virologic = kudu_virologic,
    pub_week = order_week,
    season = "2019/2020",
    prob_no_onset = filter(prob_no_onset, season == "2019/2020")
  ) 

kudu_pred$Value <- trimws(format(kudu_pred$Value, scientific = FALSE))
  
write_csv(kudu_pred,
          path = paste0("Forecasts/Live/2019-2020/Subtype Historical Average/EW", EW_paste, ".csv"))
  
write_csv(kudu_pred, 
          path = paste0(kudu_flusight_path, "/EW", EW_paste,
                        "-", substr(epiweek, 1, 4), "-Protea_Kudu.csv"))

write_csv(kudu_pred, 
          path = paste0(kudu_shiny_path, "/EW", EW_paste, ".csv"))


##### Springbok #####
best_k_cv <- readRDS("data/CV_Fourier_terms.RDS")
best_arima_cv <- readRDS("Data/CV_ARIMA_terms.RDS")
best_covar_cv <- readRDS("Data/CV_covar_terms.RDS")
gtrend_arima_fits <- readRDS('Data/gtrend_arima_fits.RDS') %>%
  filter(season == '2019/2020')

fits <- filter(flu_data_merge, year <= 2019,
               season != "2019/2020") %>%
  # Nest by location
  nest(data = c(epiweek, ILI, year, week, season, cum_h1per_6wk, cum_h3per_6wk, 
                cum_bper_6wk, region_hits, hits, order_week)) %>%
  # Create time series of ILI
  mutate(data = map(data,
                    ~ mutate(.x,
                             ILI = ts(ILI, frequency = 52, start = c(2006, 40))))) %>%
  # Merge best lambda, Fourier K value, and ARIMA structure by location
  left_join(best_k_cv, by = "location") %>%
  left_join(best_arima_cv, by = "location") %>%
  left_join(best_covar_cv, by = "location") %>%
  # Fit models
  mutate(fit = case_when(
    model == "ARIMA only" ~
      pmap(
        list(data, arima_1, arima_2, arima_3, K),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = fourier(..1$ILI, K = ..5),
                lambda = 0)
      ),
    model == "National Gtrends" ~
      pmap(
        list(data, arima_1, arima_2, arima_3, K),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$hits) %>%
                  as.matrix(),
                lambda = 0)
      ),
    model == "Regional Gtrends" ~
      pmap(
        list(data, arima_1, arima_2, arima_3, K),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$region_hits) %>%
                  as.matrix(),
                lambda = 0)
      ),
    model == "FluVirus" ~
      pmap(
        list(data, arima_1, arima_2, arima_3, K),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$cum_h1per_6wk, ..1$cum_h3per_6wk) %>%
                  as.matrix(),
                lambda = 0)
      ),
    model == "National Gtrends & FluVirus" ~
      pmap(
        list(data, arima_1, arima_2, arima_3, K),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$hits, ..1$cum_h1per_6wk,
                             ..1$cum_h3per_6wk) %>%
                  as.matrix(),
                lambda = 0)
      ),
    model == "Regional Gtrends & FluVirus" ~
      pmap(
        list(data, arima_1, arima_2, arima_3, K),
        ~ Arima(..1$ILI, order = c(..2, ..3, ..4),
                xreg = cbind(as.data.frame(fourier(..1$ILI, K = ..5)),
                             ..1$region_hits, ..1$cum_h1per_6wk,
                             ..1$cum_h3per_6wk) %>%
                  as.matrix(),
                lambda = 0)
      )
  )) %>%
  mutate(season = "2019/2020",
         prev_season = "2018/2019",
         week = EW,
         order_week = order_week,
         epiweek = epiweek) %>%
  select(-data)

# Create and save forecast files
springbok_path <- "Forecasts/Live/2019-2020/Dynamic Harmonic Model"

springbok_flusight_path <- 
  "../cdc-flusight-ensemble/model-forecasts/real-time-component-models/Protea_Springbok" 

springbok_cdc_path <- "CDC Submissions/2019-2020/Springbok"

springbok_shiny_path <- 
  '../ForecastShiny/RawData/2019-2020/NatReg/Dynamic Harmonic Model'

springbok_pred <- fits %>%
  # Join backfill densities
  left_join(filter(ili_backfill_densities, week == order_week) %>%
              select(location, season, backfill),
            by = c('location', 'season')) %>%
  # Join Google Trend Arima models
  left_join(filter(gtrend_arima_fits, location == "US National") %>%
              select(season, nat_model = arima_model) ,
            by = "season") %>%
  left_join(select(gtrend_arima_fits, location, season, region_model = arima_model),
            by = c("season", "location")) %>%
  mutate(
    pred_data = pmap(list(season, week, location, epiweek), 
                     ~ filter(flu_data_merge, year <= as.numeric(substr(..1, 6, 9)),
                              season != paste0(substr(..1, 6, 9), "/",
                                               as.numeric(substr(..1, 6, 9)) + 1),
                              season != ..1 | order_week %in% 40:..2,
                              location == ..3) %>%
                       select(epiweek, location, ILI, season, week, order_week,
                              hits, region_hits, cum_h1per_6wk, cum_h3per_6wk) %>%
                       left_join(select(ili_init_pub_list[[paste(..4)]],
                                        ILI, epiweek, location),
                                 by = c("epiweek", "location")) %>%
                       mutate(ILI = ifelse(is.na(ILI.y), ILI.x, ILI.y),
                              hits = ts(hits, frequency = 52, start = c(2006, 40)),
                              region_hits = ts(region_hits, frequency = 52, start = c(2006, 40))) %>%
                       select(-ILI.x, -ILI.y)),
    # Set up Fourier data for forecasting
    xreg = case_when(
      model == "National Gtrends" ~
        map(pred_data,
            ~ data.frame(hits = .$hits) %>%
              as.matrix()),
      model == "Regional Gtrends" ~
        map(pred_data, 
            ~ data.frame(region_hits = .$region_hits) %>%
              as.matrix()),
      model == "FluVirus" ~ 
        map(pred_data, 
            ~ data.frame(cum_h1per_6wk = .$cum_h1per_6wk,
                         cum_h3per_6wk = .$cum_h3per_6wk) %>%
              as.matrix()),
      model == "National Gtrends & FluVirus" ~ 
        map(pred_data,
            ~ data.frame(hits = .$hits,
                         cum_h1per_6wk = .$cum_h1per_6wk,
                         cum_h3per_6wk = .$cum_h3per_6wk) %>%
              as.matrix()),
      model == "Regional Gtrends & FluVirus" ~ 
        map(pred_data, 
            ~ data.frame(region_hits = .$region_hits,
                         cum_h1per_6wk = .$cum_h1per_6wk,
                         cum_h3per_6wk = .$cum_h3per_6wk) %>%
              as.matrix())
    ),
    max_week = ifelse(season == "2014/2015", 53, 52),
    # Create data frame of xreg terms for forecasting
    gtrend_forecast = pmap(
      list(pred_data, nat_model, location, season, week, max_week),
      ~ tibble(hits = c(flu_data_merge %>%
                          filter(location == ..3, season == ..4, 
                                 order_week == ..5 + 1) %>%
                          pull(hits),
                        tail(forecast(..1$hits, model = ..2,
                                      h = ..6 - 14 - nrow(..1[..1$season == ..4, ]))$mean, -1))) %>%
        as.matrix()
    ),
    reg_gtrend_forecast = pmap(
      list(pred_data, region_model, location, season, week, max_week),
      ~ tibble(region_hits = c(flu_data_merge %>%
                                 filter(location == ..3, season == ..4, 
                                        order_week == ..5 + 1) %>%
                                 pull(region_hits),
                               tail(forecast(..1$region_hits, model = ..2,
                                             h = ..6 - 14 - nrow(..1[..1$season == ..4, ]))$mean, -1))) %>%
        as.matrix()
    ),
    h1_per_forecast = pmap(
      list(pred_data, season, max_week),
      ~ tibble(cum_h1per_6wk = rep(last(..1$cum_h1per_6wk), ..3 - 14 - 
                                     nrow(..1[..1$season == ..2, ]))) %>%
        as.matrix()
    ),
    h3_per_forecast = pmap(
      list(pred_data, season, max_week),
      ~ tibble(cum_h3per_6wk = rep(last(..1$cum_h3per_6wk), ..3 - 14 - 
                                     nrow(..1[..1$season == ..2, ]))) %>%
        as.matrix()
    ),
    forecast_xreg = case_when(
      model == "National Gtrends" ~ gtrend_forecast,
      model == "Regional Gtrends" ~ reg_gtrend_forecast,
      model == "FluVirus" ~ map2(h1_per_forecast, h3_per_forecast,
                                 ~ cbind(.x, .y)),
      model == "National Gtrends & FluVirus" ~
        pmap(list(gtrend_forecast, h1_per_forecast, h3_per_forecast),
             ~ cbind(..1, ..2, ..3)),
      model == "Regional Gtrends & FluVirus" ~
        pmap(list(reg_gtrend_forecast, h1_per_forecast, h3_per_forecast),
             ~ cbind(..1, ..2, ..3))
    ),
    pred_results = pmap(
      list(fit, pred_data, location, season, K, max_week, backfill,
           xreg, forecast_xreg),
      ~ fit_to_forecast(object = ..1,
                        k = ..5,
                        pred_data = ..2,
                        location = ..3,
                        season = ..4,
                        max_week = ..6,
                        backfill = ..7,
                        xreg = ..8,
                        forecast_xreg = ..9,
                        npaths = 1000)
    )
  ) %>%
  select(location,pred_results)%>%
  unnest(col = c(pred_results)) %>%
  select(-forecast_week) %>%
  bind_rows(generate_point_forecasts(.)) %>%
  select(Location = location, Target = target, Type = type, 
         Unit = unit, Bin_start_incl = bin_start_incl, 
         Bin_end_notincl = bin_end_notincl, Value = value)

write_csv(springbok_pred, path = paste0(springbok_path, "/EW", EW_paste, ".csv"))

write_csv(springbok_pred, path = paste0(springbok_shiny_path, "/EW", EW_paste, ".csv"))

write_csv(springbok_pred, 
          path = paste0(springbok_flusight_path, "/EW", EW_paste,
                        "-", substr(epiweek, 1, 4), "-Protea_Springbok.csv"))

write_csv(springbok_pred, 
          path = paste0(springbok_cdc_path, "/EW", EW_paste,
                        "-Protea_Springbok-", Sys.Date(), ".csv"))

##### Cheetah #####

cheetah_flusight_path <- 
  "../cdc-flusight-ensemble/model-forecasts/real-time-component-models/Protea_Cheetah"

cheetah_cdc_path <- "CDC Submissions/2019-2020/Cheetah"

cheetah_shiny_path <- 
  '../ForecastShiny/RawData/2019-2020/NatReg/Ensemble'  

# List forecast files to be weighted
model_files <- list.files(path = "Forecasts/Live/2019-2020", recursive = TRUE, full.names = T)
model_files <- model_files[!grepl('ens', model_files)]

# List weighting schemes to be calculated
weight_files <- list.files(path = "Weights")
weight_types <- sub(
  pattern = "-weights.csv", 
  replacement = "",
  sub(pattern = "Weights/",
      replacement = "",
      weight_files)
)


for(j in 1:length(weight_files)) {
  stacking_weights <- read.csv(paste0("weights/", weight_files[j]), 
                               stringsAsFactors=FALSE) %>%
    filter(season == "2019/2020")
  stacked_name <- sub(pattern = ".csv", replacement = "", weight_files[j])
  
  # If weights by week included, reset to MMWR week
  if (any(grepl("Model.Week", names(stacking_weights)))) {
    stacking_weights <- mutate(
      stacking_weights,
      `Model.Week` = week_reset(`Model.Week`, season),
      ew = paste0("EW", str_pad(`Model.Week`, 2, "left", 0))
    ) %>%
      select(-`Model.Week`) 
  }
  
  wt_subset <- select(stacking_weights, -season)
    
  dir.create(file.path("Forecasts/Live/2019-2020",
                       paste0("ens-", stacked_name)), 
             showWarnings = FALSE)
  
  ## identify the "EWXX-YYYY" combos for files given season
  # first_year <- 2018
  # first_year_season_weeks <- 40:52
  # week_names <- c(paste0("EW", first_year_season_weeks),
  #                 paste0("EW", str_pad(1:20, 2, "left", pad = 0)))
    
  this_week <- paste0("EW", EW_paste)
  
  if (any(grepl("ew", names(wt_subset)))) {
    wt_sub_subset <- filter(wt_subset, ew == this_week) %>%
      select(-ew)
  } else {
    wt_sub_subset <- wt_subset
  }
  
  ## stack models, save ensemble file
  files_to_stack <- model_files[grepl(this_week, model_files)]
  
  file_df <- data.frame(
    file = files_to_stack, 
    model_id = word(files_to_stack, start = 4, end = 4, sep = "/"),
    stringsAsFactors = FALSE)
  
  stacked_entry <- stack_forecasts(file_df, wt_sub_subset) %>%
    select(location, target, type, unit, bin_start_incl, bin_end_notincl, value)
  stacked_file_name <- paste0(
    "Forecasts/Live/2019-2020/ens-",
    stacked_name, "/", this_week, ".csv"
  )
  write.csv(stacked_entry, file=stacked_file_name, 
            row.names = FALSE, quote = FALSE)
  
  # Save Cheetah model to FSN and CDC paths as well
  if(stacked_name == "month-target-based-weights") {
    write.csv(stacked_entry, 
              file = paste0(cheetah_flusight_path, "/", this_week,
                            "-", substr(epiweek, 1, 4), "-Protea_Cheetah.csv"),
              row.names = FALSE, quote = FALSE)
    
    write.csv(stacked_entry, 
              file = paste0(cheetah_cdc_path, "/", this_week,
                            "-Protea_Cheetah-", Sys.Date(), ".csv"),
              row.names = FALSE, quote = FALSE)
    
    write.csv(stacked_entry, 
              file = paste0(cheetah_shiny_path, "/EW", EW_paste, ".csv"),
              row.names = FALSE, quote = FALSE)
  }
}

# Recreate README with updated forecasts
rmarkdown::render("README.Rmd")

Sys.time()
