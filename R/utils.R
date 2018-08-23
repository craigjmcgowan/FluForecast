require(FluSight)
require(dplyr)


# Calculate scores from list of submissions ---------
calc_scores <- function(subs, truth, season, exclude = FALSE,
                        eval = FALSE, eval_period = NULL) {
  
  max_MMWR <- ifelse(season %in% 
                       c("1997/1998", "2003/2004", "2008/2009", "2014/2015"),
                     53, 52)
  
  scores <- data.frame()
  names(forecasts_1213$`Subtype Historical Average`)
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
    mutate(skill = exp(score)) %>%
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

