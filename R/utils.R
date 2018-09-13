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

read_forecasts <- function(dir, these_weeks = NULL) {
  
  teams <- list.dirs(path = dir, recursive = F)
  
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
          verify_entry(subs[[team_name]][[week]])
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


# Create evaluation period for particular year ----------
create_eval_period <- function(ILI, truth, season) {
  
  max_MMWR <- ifelse(season %in% 
                       c("1997/1998", "2003/2004", "2008/2009", "2014/2015"),
                     53, 52)
  
  # Boundaries of MMWR weeks ILINet was above baseline
  wks_abv_baseline <- ILI %>%
    left_join(FluSight::past_baselines %>%
                filter(year == as.numeric(substr(season, 1, 4))),
              by = "location") %>%
    group_by(location) %>%
    filter(ILI >= value) %>%
    summarize(end_week = last(week)) %>%
    left_join(truth %>% filter(target == "Season onset") %>% 
                mutate(start_week = case_when(
                  bin_start_incl == "none" ~ 43,
                  TRUE ~ as.numeric(bin_start_incl)
                )) %>%
                select(location, start_week),
              by = "location")
  
  # Seasonal target bounds for evaluation period
  seasonal_eval_period <- truth %>%
    # Onset weekly bins
    filter(target == "Season onset") %>%
    mutate(start_week = 43,
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
                mutate(start_week = 43,
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
                mutate(start_week = 43,
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
      bin_start_incl == "none" ~ 43,
      as.numeric(bin_start_incl) - 4 < 1 ~ as.numeric(bin_start_incl) + 48,
      as.numeric(bin_start_incl) - 4 < 43 & as.numeric(bin_start_incl) > 17 ~ 43,
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
create_eval_scores <- function(scores, bounds) {
  
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
fit_to_forecast <- function(object, xreg, pred_data, location, npaths = 1000, max_week = 52, ...){

  # Simulate output
  sim_output <- sample_predictive_trajectories_arima(
    object, 
    h = max_week - 17 - nrow(filter(pred_data, season == "2017/2018")),
    xreg = xreg,
    npaths = npaths
  )
  
  
  # Calculate forecast probabilities
  current_week <- ifelse(last(pred_data$week) < 40, last(pred_data$week) + max_week,
                         last(pred_data$week))
  season_max_ili <- filter(pred_data, season == "2017/2018") %>% pull(ILI) %>%
    max() %>% round(1)
  
  season_max_week <- filter(pred_data, season == "2017/2018") %>%
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
  calc_onset <- rbind(
    matrix(filter(pred_data, season == "2017/2018") %>%
             mutate(ILI = as.numeric(ILI)) %>% pull(ILI),
           nrow = filter(pred_data, season == "2017/2018") %>% nrow(),
           ncol = npaths),
    t(sim_output)
  ) %>% as.tibble()

  onsets <- apply(calc_onset, 2, function(x) {
    temp <- tibble(week = 40:(max_week + 22),
                   location = location,
                   ILI = x)
    try(create_onset(temp, region = location, year = 2017)$bin_start_incl, silent = TRUE)
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

  return(forecast_results)
}
