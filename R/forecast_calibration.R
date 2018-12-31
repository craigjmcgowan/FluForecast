require(tidyverse)

# Calculate calibration for weekly targets
forecast_calibration <- function(df, percent_bin, exact = TRUE,
                                 ssn_group = FALSE, chal_group = FALSE,
                                 tgt_group = FALSE, loc_group = FALSE, 
                                 team_group = FALSE, ret_teams = c(), 
                                 ens_teams = c(), repeat_teams = c(),
                                 levels = c()) {
  
  
  # Create flu/dengue levels if none specified
  if(is.null(levels)) levels = c("HHS Region 1", "HHS Region 2", 
                                 "HHS Region 3", "HHS Region 4",
                                 "HHS Region 5", "HHS Region 6",
                                 "HHS Region 7", "HHS Region 8",
                                 "HHS Region 9", "HHS Region 10",
                                 "US National", "Iquitos", "San Juan")
  
  # Remove any invalid forecasts
  df <- filter(df, !is.na(value), value <= 1, value >= 0) %>%
    # Remove unnecessary variables
    select(-bin_start_incl, -bin_end_notincl)
  
  # Check to see if all seasons have teams provided if stats requested by team type
  if (!is.null(ret_teams) && 
      (!all.equal(unique(df$season), unique(ret_teams$season)))) {
    warning("In ret_teams, not all seasons have teams identified")
  }
  if (!is.null(ens_teams) && 
      (!all.equal(unique(df$season), unique(ens_teams$season)))) {
    warning("In ens_teamss, not all seasons have teams identified")
  }
  if (!is.null(repeat_teams) && 
      length(repeat_teams) != length(unique(df$season))) {
    stop("Provide list of vectors of names of repeating teams for each year of submissions")
  }
  
  # Change location into a factor for plotting
  df <- df %>%
    mutate(location = factor(location, 
                             levels = levels))
    
  # If not exact, combine all "true" and "non-true" probabilities together
  if(exact == FALSE){
    
    df <- group_by(df, location, target, forecast_week, true)
    
    if(isTRUE(chal_group)) df <- group_by(df, challenge, add = TRUE)
    if(isTRUE(ssn_group)) df <- group_by(df, season, add = TRUE)
    if(isTRUE(tgt_group)) df <- group_by(df, target, add = TRUE)
    if(isTRUE(loc_group)) df <- group_by(df, location, add = TRUE)
    if(isTRUE(team_group)) df <- group_by(df, team, add = TRUE)
    if(!is.null(ret_teams)) df <- group_by(df, ret, add = TRUE)
    if(!is.null(ens_teams)) df <- group_by(df, ens, add = TRUE)
    
    df <- df %>%
      mutate(value = sum(value)) %>%
      unique() %>%
      ungroup()
      
  }
  
  if(!is.null(ens_teams)) {
    df <- right_join(df, ens_teams, by = c("season", "team"))
  }
  if(!is.null(ret_teams)) {
    df <- right_join(df, ret_teams, by = c("season", "team"))
  }

  # If comparing teams across seasons for repeat teams, only keep values from those teams
  if (!is.null(repeat_teams)) {
    df <- filter(df, team %in% repeat_teams)
  }
  
  # Group probabilities according to defined sensitivity of analyses (10%, 5%, etc)
  df <- df %>%
    mutate(grouping = floor(value/percent_bin),
           # Add predicted prob of exactly 1 into top group
           grouping = case_when(
             grouping == max(grouping) ~ max(grouping) - 1,
             TRUE ~ grouping
           )) %>%
    # Remove any negative predicted probabilities
    filter(grouping >= 0)
  
  # Group according to location, target, and team as needed
  df <- group_by(df, grouping)
  if(isTRUE(chal_group)) df <- group_by(df, challenge, add = TRUE)
  if(isTRUE(ssn_group)) df <- group_by(df, season, add = TRUE)
  if(isTRUE(tgt_group)) df <- group_by(df, target, add = TRUE)
  if(isTRUE(loc_group)) df <- group_by(df, location, add = TRUE)
  if(isTRUE(team_group)) df <- group_by(df, team, add = TRUE)
  if(!is.null(ret_teams)) df <- group_by(df, ret, add = TRUE)
  if(!is.null(ens_teams)) df <- group_by(df, ens, add = TRUE)

  df %>%
    summarize(percent = mean(true, na.rm = T),
              mean_pred_prob = mean(value, na.rm = T),
              n = sum(!is.na(true))) %>%
    ungroup() %>%
    mutate(min_percent = grouping * percent_bin,
           mid_percent = min_percent + percent_bin / 2,
           max_percent = min_percent + percent_bin,
           width = percent_bin) %>%
    # Calculate Jeffreys CI
    mutate(lower = qbeta(0.025, 0.5 + (n*percent), n - (n*percent) + 1/2),
           upper = qbeta(0.975, 0.5 + (n*percent), n - (n*percent) + 1/2))
  
}