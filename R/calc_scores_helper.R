require(FluSight)
require(dplyr)

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



