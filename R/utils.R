require(FluSight)
require(dplyr)

# Helper functions for reading in EpiData -------
pull_curr_epidata <- function(start, end) {
  Epidata$fluview(list('nat', 'hhs1', 'hhs2', 'hhs3', 'hhs4', 'hhs5',
                       'hhs6', 'hhs7', 'hhs8', 'hhs9', 'hhs10'),
                  list(Epidata$range(start, end)))$epidata %>%
    modify_depth(2, function(x) ifelse(is.null(x), NA, x)) %>%
    bind_rows()
}

# Helper function for reading in weekly publications of EpiData -----
pull_initpub_epidata <- function(issue) {

  start_week <- ifelse(as.numeric(substr(issue - 40, 5, 6)) < 52 & 
                         as.numeric(substr(issue - 40, 5, 6)) > 0,
                       issue - 40,
                       issue - 40 - 48) # Converts to appropriate week of previous year
  
  Epidata$fluview(list('nat', 'hhs1', 'hhs2', 'hhs3', 'hhs4', 'hhs5',
                       'hhs6', 'hhs7', 'hhs8', 'hhs9', 'hhs10'),
                  list(Epidata$range(start_week, issue)),
                  issues = list(issue))$epidata %>%
    modify_depth(2, function(x) ifelse(is.null(x), NA, x)) %>%
    bind_rows()
}


# calculate ratio of two sets of Google Trends data -------
US_flu_ratio <- function(gdata1, gdata2) {
  inner_join(gdata1 %>%
               select(date, old_hits = hits) %>%
               filter(MMWRweek(date)[[2]] > 40 | MMWRweek(date)[[2]] < 20),
             gdata2 %>%
               select(date, new_hits = hits) %>%
               filter(MMWRweek(date)[[2]] > 40 | MMWRweek(date)[[2]] < 20),
             by = "date") %>%
    mutate(ratio = old_hits / new_hits) %>%
    pull(ratio) %>%
    mean()
}

# Helper function for pulling in Google Trends data for a region -------
fetch_gtrend <- function(location, query = 'flu') {
  require(gtrendsR)
  
  # Fix location to be in gtrend format if it isn't already
  if(!location %in% c("US", state.abb, paste0("US-", state.abb)))
    stop("Invalid location provided. Must be 'US' or a two-digit state abbreviation")
  
  if(location %in% state.abb)
    location <- paste0("US-", location)
  
  flu_0407 <- gtrends(keyword = query,
                      geo = location,
                      time = paste(MMWRweek2Date(2004, 40), MMWRweek2Date(2007, 40)))$interest_over_time %>%
    mutate_at(vars(hits), as.character) %>%
    mutate(hits = case_when(hits %in% c("0", "<1") ~ '0.1',
                            TRUE ~ hits),
           hits = as.numeric(hits))
  
  flu_0611 <- gtrends(keyword = query,
                      geo = location,
                      time = paste(MMWRweek2Date(2006, 40), MMWRweek2Date(2011, 40)))$interest_over_time %>%
    mutate_at(vars(hits), as.character) %>%
    mutate(hits = case_when(hits %in% c("0", "<1") ~ '0.1',
                            TRUE ~ hits),
           hits = as.numeric(hits))
  
  flu_1015 <- gtrends(keyword = query,
                      geo = location,
                      time = paste(MMWRweek2Date(2010, 40), MMWRweek2Date(2015, 40)))$interest_over_time %>%
    mutate_at(vars(hits), as.character) %>%
    mutate(hits = case_when(hits %in% c("0", "<1") ~ '0.1',
                            TRUE ~ hits),
           hits = as.numeric(hits))
  
  flu_1419 <- gtrends(keyword = query,
                      geo = location,
                      time = paste(MMWRweek2Date(2014, 40), MMWRweek2Date(2019, 40)))$interest_over_time %>%
    mutate_at(vars(hits), as.character) %>%
    mutate(hits = case_when(hits %in% c("0", "<1") ~ '0.1',
                            TRUE ~ hits),
           hits = as.numeric(hits))
  
  flu_1820 <- gtrends(keyword = query,
                      geo = location,
                      paste(MMWRweek2Date(2018, 40), Sys.Date()))$interest_over_time %>%
    mutate_at(vars(hits), as.character) %>%
    mutate(hits = case_when(hits %in% c("0", "<1") ~ '0.1',
                            TRUE ~ hits),
           hits = as.numeric(hits))
  
  # Set up ratios to normalize all Gtrends data to scale from 2010-2015
  gratio_0607 <- US_flu_ratio(flu_0407, flu_0611)
  gratio_1011 <- US_flu_ratio(flu_0611, flu_1015)
  gratio_1415 <- US_flu_ratio(flu_1015, flu_1419)
  gratio_1819 <- US_flu_ratio(flu_1419, flu_1820)
  
  
  # Merge Gtrends data and rescale to 2010-2015 scale
  temp <- filter(flu_0407, !date %in% flu_0611$date) %>%
    mutate(hits = hits / gratio_0607) %>%
    bind_rows(flu_0611) %>%
    filter(!date %in% flu_1015$date) %>%
    mutate(hits = hits / gratio_1011) %>%
    bind_rows(flu_1015) %>%
    filter(!date %in% flu_1419$date) %>%
    mutate(hits = hits / gratio_1415) %>%
    bind_rows(flu_1419) %>%
    filter(!date %in% flu_1820$date) %>%
    mutate(hits = hits / gratio_1819) %>%
    bind_rows(flu_1820) %>%
    select(date, hits) %>%
    do({tz(.$date) <- "America/New_York"; .}) %>%
    mutate(week = MMWRweek(date)[[2]],
           year = MMWRweek(date)[[1]],
           season = ifelse(week >= 40,
                           paste0(year, "/", year + 1),
                           paste0(year - 1, "/", year)),
           hits = case_when(near(hits, 0) ~ 1,
                            TRUE ~ hits))
}

# Functions to keep weeks in order and reset them --------
week_inorder <- function(week, season){
  case_when(as.numeric(week) < 40 & 
              season %in% c("1997/1998", "2003/2004", "2008/2009", "2014/2015") ~ 
              as.numeric(week) + 53,
            as.numeric(week) < 40 ~ as.numeric(week) + 52,
            TRUE ~ as.numeric(week))
}

week_reset <- function(week, season) {
  case_when(as.numeric(week) > 53 & 
              season %in% c("1997/1998", "2003/2004", "2008/2009", "2014/2015") ~ 
              as.numeric(week) - 53,
            as.numeric(week) > 52 & 
              !season %in% c("1997/1998", "2003/2004", "2008/2009", "2014/2015") ~
              as.numeric(week) - 52,
            TRUE ~ as.numeric(week))
}


# Create pseudo-onsets from calculated baselines --------
create_pseudo_onset <- function(weekILI) {
  
  # Save maximum MMWR week in season being analyzed
  maxMMWR <- max(weekILI$week)
  
  # Add 52/53 to weeks in new year to keep weeks in order
  weekILI$week[weekILI$week < 40] <-
    as.integer(weekILI$week[weekILI$week < 40] + maxMMWR)
  
  # Check to see if 3 weeks above baseline have passed
  j <- 0  # Counter for weeks above peak
  for (i in head(weekILI$week, n = 1):tail(weekILI$week, n = 1)) {
    if (weekILI$ILI[weekILI$week == i] >= weekILI$baseline[weekILI$week == i]) {
      j <- j + 1
    } else {
      j <- 0
    }
    if (j == 3) {
      onset <- i - 2
      break
    }
    if (i == tail(weekILI$week, n = 1)) {
      onset <- "none"
    }
  }
  
  # If onset week > 52, reset to MMWR week
  if (is.numeric(onset) && onset > maxMMWR) {
    onset <- onset - maxMMWR
  }
  
  onset_truth <- data.frame(target = "Season onset",
                            location = weekILI$location[1],
                            forecast_week = as.integer(NA),
                            bin_start_incl = as.character(onset),
                            stringsAsFactors = FALSE) %>%
    mutate(bin_start_incl = trimws(replace(bin_start_incl,
                                           !is.na(bin_start_incl) & bin_start_incl != "none",
                                           format(round(as.numeric(
                                             bin_start_incl[!is.na(bin_start_incl) & bin_start_incl != "none"])
                                             , 1), nsmall = 1, trim = T))))
  
  
  return(onset_truth)
}  


# Functions to read forecasts from directory into a list --------

read_forecasts <- function(dir, challenge = 'ilinet', these_weeks = c(),
                           these_teams = c()) {
  require(stringr)
  
  teams <- list.dirs(path = dir, recursive = F)
  
  # Only read forecasts for certain teams if called
  if(!is.null(these_teams))
    teams <- teams[str_extract(teams, "[^/]+$") %in% these_teams]
  
  subs <- list()
  
  for (this_team in teams) {
    
    # Extract team name from folder
    team_name <- str_extract(this_team, "[^/]+$")
    # Vector of submission files
    files <- list.files(path = this_team, pattern = "*.csv")
    # Only keep submission files from weeks specified
    if (!is.null(these_weeks)) {
      files <- files[str_extract(files, "EW[0-9]{2}") %in% these_weeks]
    }
    
    for (this_file in files) {
      
      week <- str_extract(this_file, "EW[0-9]{1,2}")
      
      tryCatch(
        {
          subs[[team_name]][[week]] <- read_entry(paste0(this_team,"/", this_file))
        },
        error = function(cond) {
          message(paste("Errors reading in submission for", team_name))
          message(cond)
          return(NA)
        },
        warning = function(cond) {
          message(paste("Warnings reading in submission for", team_name))
          message(cond)
          return(NA)
        },
        finally = message(paste("Submission read in for", team_name))
      )
      
      tryCatch(
        {
          verify_entry(subs[[team_name]][[week]], challenge = challenge)
        },
        error = function(cond) {
          message(paste("Errors in verifying submission for", team_name))
          message(cond)
          return(NA)
        },
        warning = function(cond) {
          message(paste("Warnings in verifying submission for", team_name))
          message(cond)
          return(NA)
        },
        finally = message(paste(week, "submission processed for", team_name))
      )
      
      # Remove invalid probabilities and normalize other probabilities.
      tryCatch(
        {
          subs[[team_name]][[week]] <- 
            remove_invalid(subs[[team_name]][[week]])
          
          subs[[team_name]][[week]] <- 
            normalize_probs(subs[[team_name]][[week]], ignore_invalid = TRUE)
        },
        error = function(cond) {
          message(paste("Errors in normalizing submission for", team_name))
          message(cond)
          return(NA)
        },
        warning = function(cond) {
          message(paste("Warnings in normalizing submission for", team_name))
          message(cond)
          return(NA)
        },
        finally = message(paste(week, "Probabilities normalized for", team_name))
      )
    }
  }
  
  
  return(subs)
}

subs_by_week <- function(subs) {
  outlist <- list()
  
  for (this_week in names(subs[[1]])){
    week_data <- data.frame()
    for (this_team in names(subs)) {
      # Check if entry for that week is missing
      if (!is.null(subs[[this_team]][[this_week]])) { 
        week_data <- subs[[this_team]][[this_week]] %>%
          select(-forecast_week) %>%
          mutate(team = this_team) %>%
          bind_rows(week_data)
      }
    }
    outlist[[this_week]] <- week_data %>%
      arrange(team)
  }
  
  return(outlist)
}

# Calculate scores ------
calc_scores <- function(subs, truth, season, exclude = FALSE,
                        eval = FALSE, eval_period = NULL) {
  
  max_MMWR <- ifelse(season %in% 
                       c("1997/1998", "2003/2004", "2008/2009", "2014/2015"),
                     53, 52)
  
  scores <- data.frame()
  
  for (this_team in names(subs)) {
    
    # Determine which weeks to score
    # exclude = F will count missing weeks as -10
    # exclude = T will only score forecasts that teams have submitted
    if (exclude == FALSE) {
      weeks <- c("EW01", "EW02", "EW03", "EW04", "EW05", "EW06", "EW07", "EW08",
                 "EW09", "EW10", "EW11", "EW12", "EW13", "EW14", "EW15", "EW16",
                 "EW17", "EW18", "EW19", "EW20", "EW40", "EW41", "EW42", "EW43",
                 "EW44", "EW45", "EW46", "EW47", "EW48", "EW49", "EW50", "EW51",
                 "EW52")
    } else weeks <- names(subs[[this_team]])
    for (this_week in weeks) { 
      # Check for missing forecast - assign -10 if missing
      if (is.null(subs[[this_team]][[this_week]])) {
        these_scores <- expand.grid(location = c("US National", "HHS Region 1", 
                                                 "HHS Region 2", "HHS Region 3",
                                                 "HHS Region 4", "HHS Region 5", 
                                                 "HHS Region 6", "HHS Region 7",
                                                 "HHS Region 8", "HHS Region 9", 
                                                 "HHS Region 10"),
                                    target = c("Season onset", "Season peak percentage",
                                               "Season peak week", "1 wk ahead",
                                               "2 wk ahead", "3 wk ahead", 
                                               "4 wk ahead"),
                                    stringsAsFactors = FALSE)
        these_scores$score <- -10
        these_scores$forecast_week <- as.numeric(gsub("EW", "", this_week))
        these_scores$team <- this_team
      } else {
        # Score receieved entry
        these_scores <- score_entry(subs[[this_team]][[this_week]], truth)
        these_scores$team <- this_team
      }
      # Bind entry scores together
      scores <- bind_rows(scores, these_scores)
    }
  }
  
  # Create forecast skill metric and sort
  scores <- scores %>% 
    mutate(season = season) %>%
    arrange(team, forecast_week)
  
  # Return scores for evaluation period if requested
  if (isTRUE(eval)) {
    scores <- scores %>%
      left_join(eval_period, by = c("location", "target")) %>%
      # Get all weeks in order - deal w/ New Year transition
      mutate(forecast_week_order = ifelse(forecast_week < 40, 
                                          forecast_week + max_MMWR, forecast_week),
             start_week_order = ifelse(start_week < 40, 
                                       start_week + max_MMWR, start_week),
             end_week_order = ifelse(end_week < 40, 
                                     end_week + max_MMWR, end_week)) %>%
      filter(forecast_week_order >= start_week_order &
               forecast_week_order <= end_week_order) %>%
      select(-start_week, -end_week, -forecast_week_order, 
             -start_week_order, -end_week_order)
  }
  
  return(scores)
}

# Create evaluation period for particular year ----------
create_eval_period <- function(ILI, truth, season) {
  
  max_MMWR <- ifelse(season %in% 
                       c("1997/1998", "2003/2004", "2008/2009", "2014/2015"),
                     53, 52)
  
  # Boundaries of MMWR weeks ILINet was above baseline
  wks_abv_baseline <- ILI %>%
    mutate(ILI = round(ILI, 1)) %>%
    left_join(FluSight::past_baselines %>%
                filter(year == as.numeric(substr(season, 1, 4))),
              by = "location") %>%
    filter(ILI >= value) %>%
    group_by(location) %>%
    summarize(end_week = last(week)) %>%
    left_join(truth %>% filter(target == "Season onset") %>% 
                mutate(start_week = case_when(
                  bin_start_incl == "none" ~ 42,
                  TRUE ~ as.numeric(bin_start_incl)
                )) %>%
                select(location, start_week),
              by = "location")
  
  # Seasonal target bounds for evaluation period
  seasonal_eval_period <- truth %>%
    # Onset weekly bins
    filter(target == "Season onset") %>%
    mutate(start_week = 42,
           end_week = case_when(
             bin_start_incl == "none" ~ 18,
             TRUE ~ as.numeric(bin_start_incl) + 6
           ),
           end_week = case_when(
             end_week > max_MMWR ~ end_week - max_MMWR,
             TRUE ~ end_week
           ),
           end_week_order = case_when(
             end_week < 40 ~ end_week + max_MMWR,
             TRUE ~ end_week
           )) %>%
    # Peak week and peak percent
    bind_rows(truth %>%
                filter(target == "Season peak week") %>%
                left_join(wks_abv_baseline, by = "location") %>%
                mutate(start_week = 42,
                       end_week = case_when(
                         is.na(end_week) ~ 18,
                         TRUE ~ as.numeric(end_week) + 1
                       ),
                       end_week = case_when(
                         end_week > 18 & end_week < 40 ~ 18,
                         TRUE ~ end_week
                       ),
                       end_week_order = case_when(
                         end_week < 40 ~ end_week + max_MMWR,
                         TRUE ~ end_week
                       ))) %>%
    bind_rows(truth %>%
                filter(target == "Season peak week") %>%
                left_join(wks_abv_baseline, by = "location") %>%
                mutate(start_week = 42,
                       end_week = case_when(
                         is.na(end_week) ~ 18,
                         TRUE ~ as.numeric(end_week) + 1
                       ),
                       end_week = case_when(
                         end_week > 18 & end_week < 40 ~ 18,
                         TRUE ~ end_week
                       ),
                       end_week_order = case_when(
                         end_week < 40 ~ end_week + max_MMWR,
                         TRUE ~ end_week
                       ),
                       target = "Season peak percentage")) %>%
    select(target, location, start_week, end_week)
  
  # Week eval period 
  single_week_eval_period <- truth %>%
    filter(target == "Season onset") %>%
    mutate(start_week = case_when(
      bin_start_incl == "none" ~ 42,
      as.numeric(bin_start_incl) - 4 < 1 ~ as.numeric(bin_start_incl) + 48,
      as.numeric(bin_start_incl) - 4 < 42 & as.numeric(bin_start_incl) > 17 ~ 42,
      TRUE ~ as.numeric(bin_start_incl) - 4
    )) %>%
    left_join(wks_abv_baseline %>% select(-start_week), by = "location") %>%
    mutate(end_week = case_when(
      is.na(end_week) ~ 18,
      end_week + 3 > 18 & end_week + 3 < 40 ~ 18,
      end_week + 3 > max_MMWR ~ end_week - (max_MMWR - 3),
      TRUE ~ end_week + 3
    )) %>%
    select(-forecast_week, -bin_start_incl)
  
  week_eval_period <- bind_rows(
    single_week_eval_period %>%
      mutate(target = "1 wk ahead"),
    single_week_eval_period %>%
      mutate(target = "2 wk ahead"),
    single_week_eval_period %>%
      mutate(target = "3 wk ahead"),
    single_week_eval_period %>%
      mutate(target = "4 wk ahead")
  )
  
  # Combine eval bounds into single dataframe
  bind_rows(seasonal_eval_period, week_eval_period) %>%
    mutate(start_week = ifelse(start_week < 40, start_week + max_MMWR, start_week),
           end_week = ifelse(end_week < 40, end_week + max_MMWR, end_week)) %>%
    unique()
}


# Restrict score to evaluation period ------
create_eval_scores <- function(scores, bounds, season, ...) {
  
  max_MMWR <- ifelse(season %in% 
                       c("1997/1998", "2003/2004", "2008/2009", "2014/2015"),
                     53, 52)
  
  scores %>%
    left_join(bounds, by = c("location", "target")) %>%
    # Get all weeks in order - deal w/ New Year transition
    mutate(forecast_week_order = ifelse(forecast_week < 40, 
                                        forecast_week + max_MMWR, forecast_week),
           start_week_order = ifelse(start_week < 40, 
                                     start_week + max_MMWR, start_week),
           end_week_order = ifelse(end_week < 40, 
                                   end_week + max_MMWR, end_week)) %>%
    filter(forecast_week_order >= start_week_order &
             forecast_week_order <= end_week_order) %>%
    select(-start_week, -end_week, -forecast_week_order, 
           -start_week_order, -end_week_order)
}

# Simulate predicted trajectories from ARIMA models ------
sample_predictive_trajectories_arima <- function (object,
                                                  h = ifelse(object$arma[5] > 1, 2 * object$arma[5]),
                                                  xreg = NULL,
                                                  lambda = object$lambda,
                                                  npaths = 5000,
                                                  ...) {
  sim <- matrix(NA, nrow = npaths, ncol = h)

  for (i in 1:npaths) {
    try(sim[i, ] <- forecast:::simulate.Arima(object, 
                                              xreg = xreg,
                                              nsim = h,
                                              bootstrap = TRUE), silent = TRUE)
  }
  
  sim <- round(sim, 1)
  
  return(sim)
  
}

# Turn fitted forecast model into a CDC-style forecast ------
fit_to_forecast <- function(object, k, pred_data, location, season, week,
                            xreg = c(), forecast_xreg = c(),
                            backfill, npaths = 1000, max_week = 52, 
                            challenge = "ilinet", ...){

  # Create empty matrices to fill in while simulating
  sim_season <- matrix(nrow = npaths,
                       ncol = nrow(pred_data[pred_data$season == season, ]))
  
  sim_output <- matrix(nrow = npaths, 
                       ncol = max_week - 14 - nrow(pred_data[pred_data$season == season, ]))
  
  for(i in 1:npaths) {
  
    # Simulate backfill and create ILI time-series
    pred_data_temp <- pred_data %>%
      left_join(
        backfill %>%
          mutate(adjust = map2_dbl(data, density_bw,
                               ~ sample(.x$backfill, 1) + rnorm(1, 0, .y)),
                 season = season) %>%
          select(season, order_week = measure_week, adjust),
        by = c('season', 'order_week')
      ) %>%
      mutate(adjust = case_when(is.na(adjust) ~ 0,
                                TRUE ~ adjust),
             ILI = ILI + adjust,
             ILI = ts(ILI, frequency = 52, start = c(2006, 40)))
    
    xreg_temp = cbind(fourier(pred_data_temp$ILI, K = k), xreg)
    
    forecast_xreg_temp = cbind(fourier(pred_data_temp$ILI, K = k,
                                       h = max_week - 14 - 
                                         nrow(pred_data[pred_data$season == season, ])),
                               forecast_xreg)
    
    # Save matrix of current season estimates for use in onset later
    sim_season[i, ] <- pred_data_temp[pred_data_temp$season == season, ] %>%
                                mutate(ILI = as.numeric(ILI)) %>% pull(ILI)
    
    
    # Simulate output
    sim_output[i, ] <- sample_predictive_trajectories_arima(
        Arima(pred_data_temp$ILI, xreg = xreg_temp, model = object), 
        h = max_week - 14 - nrow(pred_data_temp[pred_data_temp$season == season, ]),
        xreg = forecast_xreg_temp,
        npaths = 1
      )
    
  }
  
  
  # Calculate forecast probabilities
  current_week <- ifelse(last(pred_data$week) < 40, last(pred_data$week) + max_week,
                         last(pred_data$week))
  season_max_ili <- pred_data[pred_data$season == season, ] %>% pull(ILI) %>%
    max() %>% round(1)
  
  season_max_week <- pred_data[pred_data$season == season, ] %>%
    mutate(ILI = as.numeric(ILI)) %>%
    filter(ILI == max(ILI)) %>%
    mutate(week = ifelse(week < 40, week + max_week, week)) %>%
    pull(week)
  
  forecast_results <- tibble()

  # Week ahead forecasts
  for (i in 1:4) {
    for (j in seq(0, 12.9, 0.1)) {
      forecast_results <- bind_rows(
        forecast_results,
        tibble(target = paste(i, "wk ahead"),
               bin_start_incl = trimws(format(round(j, 1), nsmall = 1)),
               value = sum(near(sim_output[, i],j), na.rm = T) /
                 length(na.omit(sim_output[, i])))
      )
    }
    forecast_results <- bind_rows(
      forecast_results,
      tibble(target = paste(i, "wk ahead"),
             bin_start_incl = "13.0",
             value = sum((sim_output[, i] >= 13), na.rm = T) /
               length(na.omit(sim_output[, i])))
    )
  }

  # Peak percentage forecasts
  max_ili <- apply(sim_output, 1, max, na.rm = T)

  # If max predicted is less than maximum observed so far this season, replace with that
  max_ili <- ifelse(max_ili < season_max_ili, season_max_ili, max_ili)


  for (j in seq(0, 12.9, 0.1)) {
    forecast_results <- bind_rows(
      forecast_results,
      tibble(target = "Season peak percentage",
             bin_start_incl = trimws(format(round(j, 1), nsmall = 1)),
             value = sum(near(max_ili, j))/length(max_ili))
    )
  }
  forecast_results <- bind_rows(
    forecast_results,
    tibble(target = "Season peak percentage",
           bin_start_incl = "13.0",
           value = sum((max_ili >= 13))/length(max_ili))
  )

  # Season peak week
  abv_max <- sim_output[apply(sim_output, 1, max, na.rm = T) > season_max_ili, , drop = FALSE]
  below_max <- sim_output[apply(sim_output, 1, max, na.rm = T) < season_max_ili, , drop = FALSE]
  at_max <- sim_output[near(apply(sim_output, 1, max, na.rm = T), season_max_ili), , drop = FALSE]

  peaks <- c(
    unlist(apply(abv_max,  1, function(x) which(x == max(x, na.rm = TRUE)) + current_week)),
    unlist(apply(at_max,  1, function(x) which(x == max(x, na.rm = TRUE)) + current_week)),
    rep(season_max_week, nrow(at_max) + nrow(below_max))
  )

  for (j in 40:(max_week + 20)) {
    forecast_results <- bind_rows(
      forecast_results,
      tibble(target = "Season peak week",
             bin_start_incl = trimws(format(round(j, 1), nsmall = 1)),
             value = sum(near(peaks, j))/length(peaks))
    )
  }

  # Onset week
  if(challenge == "ilinet") {
    calc_onset <- cbind(sim_season, sim_output) %>% 
      t() %>%
      as_tibble(.name_repair = 'minimal')
  
    onsets <- apply(calc_onset, 2, function(x) {
      temp <- tibble(week = 40:(max_week + 25),
                     location = location,
                     ILI = x)
      try(create_onset(temp, region = location, 
                       year = as.numeric(substr(season, 1, 4)))$bin_start_incl, silent = TRUE)
    }) %>%
      trimws()
  
    for (j in 40:(max_week + 20)) {
      forecast_results <- bind_rows(
        forecast_results,
        tibble(target = "Season onset",
               bin_start_incl = trimws(format(round(j, 1), nsmall = 1)),
               value = sum(onsets == trimws(format(round(j, 1), nsmall = 1)))/length(onsets))
      )
    }
    forecast_results <- bind_rows(
      forecast_results,
      tibble(target = "Season onset",
             bin_start_incl = "none",
             value = sum((onsets == "none"))/length(onsets))
    )
  }

  # Format in CDC style
  forecast_results <- forecast_results %>%
    mutate(type = "Bin",
           unit = case_when(
             target %in% c("Season onset", "Season peak week") ~ "week",
             TRUE ~ "percent"
             ),
           bin_start_incl = case_when(
             target %in% c("Season onset", "Season peak week") & 
               as.numeric(bin_start_incl) > max_week ~ 
               format(round(as.numeric(bin_start_incl) - max_week, 1), nsmall = 1),
             TRUE ~ bin_start_incl
           ),
           bin_end_notincl = case_when(
             bin_start_incl == "13.0" ~ "100.0",
             bin_start_incl == "none" ~ "none",
             target %in% c("Season onset", "Season peak week") ~
               trimws(format(round(as.numeric(bin_start_incl) + 1, 1),
                      nsmall = 1)),
             TRUE ~ trimws(format(round(as.numeric(bin_start_incl) + 0.1, 1),
                           nsmall = 1))
             ),
           forecast_week = last(pred_data$week),
           bin_start_incl = trimws(bin_start_incl),
           bin_end_notincl = trimws(bin_end_notincl)
    )

  # Normalize probabilities so all targets sum to 1
  forecast_results <- forecast_results %>%
    group_by(target) %>%
    mutate(value = value / sum(value)) %>%
    ungroup()
  
  return(forecast_results)
}


# Stack forecasts for ensembles -----
stack_forecasts <- function(files, stacking_weights, challenge = "ilinet") {
  require(dplyr)
  require(FluSight)
  nfiles <- nrow(files)
  
  ## retrieve expected model names from stacking_weights matrix
  model_names <- unique(stacking_weights$component_model_id)
  if(nfiles != length(model_names))
    stop("number of model_ids in weight matrix does not equal number of files")
  
  ## check that model weights sum to 1 for fixed target/location
  if("target" %in% colnames(stacking_weights)) {
    if("location" %in% colnames(stacking_weights)) {
      weight_sums <- stacking_weights %>% group_by(target, location) %>% 
        summarize(sum_weight=sum(weight)) %>% .$sum_weight
    } else {
      weight_sums <- stacking_weights %>% group_by(target) %>% 
        summarize(sum_weight=sum(weight)) %>% .$sum_weight
    }
  } else if("location" %in% colnames(stacking_weights)) {
    weight_sums <- stacking_weights %>% group_by(location) %>% 
      summarize(sum_weight=sum(weight)) %>% .$sum_weight
  } else {
    weight_sums <- stacking_weights %>% 
      summarize(sum_weight=sum(weight)) %>% .$sum_weight
  }
  if(!isTRUE(all.equal(weight_sums, rep(1, length(weight_sums)))))
    stop("weights don't sum to 1.")
  
  ## check that files are entries
  entries <- vector("list", nfiles)
  for(i in 1:nfiles) {
    entries[[i]] <- read_entry(files$file[i]) 
    verify_entry(entries[[i]], challenge = challenge)
  }
  
  ## stack distributions
  for(i in 1:length(entries)){
    entries[[i]] <- entries[[i]] %>%
      filter(type == "Bin") %>%
      ## add column with component_model_id
      mutate(component_model_id = files$model_id[i]) %>%
      ## add weights column
      left_join(stacking_weights) 
    ## rename value column
    new_value_name <- paste0(files$model_id[i], "_value")
    entries[[i]][new_value_name] <- with(entries[[i]], value)
    ## rename weight column
    new_wt_name <- paste0(files$model_id[i], "_weight")
    entries[[i]][new_wt_name] <- with(entries[[i]], weight)
    ## add weighted value column
    new_wtvalue_name <- paste0(files$model_id[i], "_weighted_value")
    entries[[i]][new_wtvalue_name] <- with(entries[[i]], weight * value)
  }
  ## drop unneeded columns
  unneeded_columns <- c("component_model_id", "value", "weight")
  slim_entries <- lapply(entries, function(x) x[!(names(x) %in% unneeded_columns)])
  ensemble_entry <- Reduce(
    f = left_join, 
    x = slim_entries) %>% 
    as_tibble() %>%
    mutate(value = rowSums(.[grep("weighted_value", names(.))], na.rm = TRUE)) %>%
    select(-contains("_value")) %>%
    select(-contains("_weight")) %>%
    select(-c(forecast_week))
  ## correct point estimates because averages break for weekly means
  corrected_point_ests <- generate_point_forecasts(ensemble_entry)
  ensemble_entry <- bind_rows(corrected_point_ests, ensemble_entry)
  return(ensemble_entry)
}

