require(tidyverse)
require(FluSight)

# Make historical densities based on provided data -----
create_historical_densities <- function(ili_df, pseudo_onsets) {
  
  # Only keep data from 2000/2001 season on, except 2009/2010
  train_ili <- filter(ili_df, season != "2009/2010",
                      year >= 2000, season != "1999/2000") %>%
    select(location, season, week, year, ILI) %>%
    mutate(order_week = week_inorder(week, season))
  
  # Generate past onsets from years with baseline data
  train_onsets <- train_ili %>%
    filter(year >= 2007, !season %in% c("2006/2007", "2009/2010"),
           (week >= 40 | week <= 20)) %>%
    group_by(location, season) %>%
    do(create_onset(., region = .$location[1], 
                    year = as.numeric(substr(.$season[1], 1, 4)))) %>%
    ungroup() %>%
    # Bind to pseudo-onsets
    bind_rows(pseudo_onsets)
  
  # Generate peaks from training data except 2008/2009 and 2009/2010 seasons
  train_peak <- train_ili %>%
    filter(season != "2008/2009", week >= 40 | week <= 20) %>%
    group_by(location, season) %>%
    do(create_peak(., location = .$location[1])) %>%
    ungroup()
  
  # Create empty object to store densities
  densities <- list()

  for(this_location in unique(train_ili$location)) {
    
    # Onset densities
    temp_onset <- train_onsets %>% 
      filter(location == this_location, bin_start_incl != "none") %>%
      mutate(bin_start_incl = case_when(
        as.numeric(bin_start_incl) < 40 & 
          season %in% c("1997/1998", "2003/2004", "2008/2009", "2014/2015") ~
          as.numeric(bin_start_incl) + 53,
        as.numeric(bin_start_incl) < 40 ~ 
          as.numeric(bin_start_incl) + 52,
        TRUE ~ as.numeric(bin_start_incl)))
    
    densities[["Season onset"]][[this_location]] <- 
      density(temp_onset$bin_start_incl, bw = "SJ")
    
    # Peak week densities
    temp_peak_wk <- train_peak %>%
      filter(location == this_location, target == "Season peak week") %>%
      mutate(bin_start_incl = case_when(
        as.numeric(bin_start_incl) < 40 & 
          season %in% c("1997/1998", "2003/2004", "2008/2009", "2014/2015") ~
          as.numeric(bin_start_incl) + 53,
        as.numeric(bin_start_incl) < 40 ~ 
          as.numeric(bin_start_incl) + 52,
        TRUE ~ as.numeric(bin_start_incl))) 
    
    densities[["Season peak week"]][[this_location]] <- 
      density(temp_peak_wk$bin_start_incl, bw = "SJ")
    
    # Peak percentage densities
    temp_peak_per <- train_peak %>%
      mutate(bin_start_incl = as.numeric(bin_start_incl)) %>%
      filter(location == this_location, target == "Season peak percentage")
    
    densities[["Season peak percentage"]][[this_location]] <- 
      density(temp_peak_per$bin_start_incl, bw = "SJ")
    
    
    # Densities for each week of the season
    for(this_week in 40:77) {
      temp_ili <- filter(train_ili, location == this_location,
                         order_week == this_week)
      
      densities[[paste(this_week)]][[this_location]] <-
        density(temp_ili$ILI, bw = "SJ")
    }
  }
  
  return(densities)
  
}

# Create unweighted historical forecasts -------
create_historical_forecast <- function(functions, pub_week, 
                                    season, prob_no_onset) {

  pred <- tibble()

  for(this_location in names(functions[["Season onset"]])) {
    
    # Seasonal targets -----
    onset <- tibble()
    pkwk <- tibble()
    
    for (i in 40:73) {
      
      # Season onset
      onset <- tibble(
        location = this_location,
        target = "Season onset",
        type = "Bin",
        unit = "week",
        bin_start_incl = paste(i),
        bin_end_notincl = paste(i + 1),
        value = integrate(functions[["Season onset"]][[this_location]],
                           i - 0.5, i + 0.5)$value
      ) %>% bind_rows(onset, .)
      
      # Season peak week
      pkwk <- tibble(
        location = this_location,
        target = "Season peak week",
        type = "Bin",
        unit = "week",
        bin_start_incl = paste(i),
        bin_end_notincl = paste(i + 1),
        value = integrate(functions[["Season peak week"]][[this_location]],
                           i - 0.5, i + 0.5)$value
      ) %>% bind_rows(pkwk, .)
      
    }
    
    # Remove extra week if not a 53 week season and reset numbers
    onset <- onset %>%
      mutate(bin_start_incl = as.character(week_reset(bin_start_incl, season)),
             bin_end_notincl = as.character(as.numeric(bin_start_incl) + 1)) %>%
    # Remove week 21 for seasons with 52 weeks
    filter(bin_start_incl != "21")
    
    pkwk <- pkwk %>%
      mutate(bin_start_incl = as.character(week_reset(bin_start_incl, season)),
             bin_end_notincl = as.character(as.numeric(bin_start_incl) + 1)) %>%
      # Remove week 21 for seasons with 52 weeks
      filter(bin_start_incl != "21")
    
    # Adjust onset probabilities to reflect possibility of no onset
    prob_no <- prob_no_onset$prob_no_onset[prob_no_onset$location == this_location]
    
    onset <- mutate(onset, value = value * (1 - prob_no)) %>%
      bind_rows(
        tibble(
          location = this_location,
          target = "Season onset",
          type = "Bin",
          unit = "week",
          bin_start_incl = "none",
          bin_end_notincl = "none",
          value = prob_no
        )
      )
    
    # ILI targets
    
    # Lower bounds
    pkper <- tibble(
      location = this_location,
      target = "Season peak percentage",
      type = "Bin",
      unit = "percent",
      bin_start_incl = paste(0),
      bin_end_notincl = paste(0.1),
      value = integrate(functions[["Season peak percentage"]][[this_location]],
                         0, 0.05)$value
    )
    
    wk1 <- tibble(
      location = this_location,
      target = "1 wk ahead",
      type = "Bin",
      unit = "percent",
      bin_start_incl = paste(0),
      bin_end_notincl = paste(0.1),
      value = integrate(functions[[paste(pub_week + 1)]][[this_location]],
                         0, 0.05)$value
    )
    
    wk2 <- tibble(
      location = this_location,
      target = "2 wk ahead",
      type = "Bin",
      unit = "percent",
      bin_start_incl = paste(0),
      bin_end_notincl = paste(0.1),
      value = integrate(functions[[paste(pub_week + 2)]][[this_location]],
                         0, 0.05)$value
    )
    
    wk3 <- tibble(
      location = this_location,
      target = "3 wk ahead",
      type = "Bin",
      unit = "percent",
      bin_start_incl = paste(0),
      bin_end_notincl = paste(0.1),
      value = integrate(functions[[paste(pub_week + 3)]][[this_location]],
                         0, 0.05)$value
    )
    
    wk4 <- tibble(
      location = this_location,
      target = "4 wk ahead",
      type = "Bin",
      unit = "percent",
      bin_start_incl = paste(0),
      bin_end_notincl = paste(0.1),
      value = integrate(functions[[paste(pub_week + 4)]][[this_location]],
                         0, 0.05)$value
    )
    
    
    for (i in seq(0.1, 12.9, 0.1)) {
      
      pkper <- tibble(
        location = this_location,
        target = "Season peak percentage",
        type = "Bin",
        unit = "percent",
        bin_start_incl = paste(i),
        bin_end_notincl = paste(i + 0.1),
        value = integrate(functions[["Season peak percentage"]][[this_location]],
                           i - 0.05, i + 0.05)$value
      ) %>% bind_rows(pkper, .)
      
      wk1 <- tibble(
        location = this_location,
        target = "1 wk ahead",
        type = "Bin",
        unit = "percent",
        bin_start_incl = paste(i),
        bin_end_notincl = paste(i + 0.1),
        value = integrate(functions[[paste(pub_week + 1)]][[this_location]],
                           i - 0.05, i + 0.05)$value
      ) %>% bind_rows(wk1, .)
      
      wk2 <- tibble(
        location = this_location,
        target = "2 wk ahead",
        type = "Bin",
        unit = "percent",
        bin_start_incl = paste(i),
        bin_end_notincl = paste(i + 0.1),
        value = integrate(functions[[paste(pub_week + 2)]][[this_location]],
                           i - 0.05, i + 0.05)$value 
      ) %>% bind_rows(wk2, .)
      
      wk3 <- tibble(
        location = this_location,
        target = "3 wk ahead",
        type = "Bin",
        unit = "percent",
        bin_start_incl = paste(i),
        bin_end_notincl = paste(i + 0.1),
        value = integrate(functions[[paste(pub_week + 3)]][[this_location]],
                           i - 0.05, i + 0.05)$value
      ) %>% bind_rows(wk3, .)
      
      wk4 <- tibble(
        location = this_location,
        target = "4 wk ahead",
        type = "Bin",
        unit = "percent",
        bin_start_incl = paste(i),
        bin_end_notincl = paste(i + 0.1),
        value = integrate(functions[[paste(pub_week + 4)]][[this_location]],
                           i - 0.05, i + 0.05)$value
      ) %>% bind_rows(wk4, .)
      
    }
    
    # Upper bounds
    pkper <- tibble(
      location = this_location,
      target = "Season peak percentage",
      type = "Bin",
      unit = "percent",
      bin_start_incl = paste(13),
      bin_end_notincl = paste(100),
      value = integrate(functions[["Season peak percentage"]][[this_location]],
                         12.95, 15)$value 
    ) %>% 
      mutate(value = ifelse(all(tail(pkper$value) == tail(pkper$value, n = 1)),
                            tail(pkper$value, n = 1), value)) %>%
      bind_rows(pkper, .)
    
    wk1 <- tibble(
      location = this_location,
      target = "1 wk ahead",
      type = "Bin",
      unit = "percent",
      bin_start_incl = paste(13),
      bin_end_notincl = paste(100),
      value = integrate(functions[[paste(pub_week + 1)]][[this_location]],
                         12.95, 15)$value
    ) %>% 
      mutate(value = ifelse(all(tail(wk1$value) == tail(wk1$value, n = 1)),
                            tail(wk1$value, n = 1), value)) %>%
      bind_rows(wk1, .)
    
    wk2 <- tibble(
      location = this_location,
      target = "2 wk ahead",
      type = "Bin",
      unit = "percent",
      bin_start_incl = paste(13),
      bin_end_notincl = paste(100),
      value = integrate(functions[[paste(pub_week + 2)]][[this_location]],
                         12.95, 15)$value
    ) %>% 
      mutate(value = ifelse(all(tail(wk2$value) == tail(wk2$value, n = 1)),
                            tail(wk2$value, n = 1), value)) %>%
      bind_rows(wk2, .)
    
    wk3 <- tibble(
      location = this_location,
      target = "3 wk ahead",
      type = "Bin",
      unit = "percent",
      bin_start_incl = paste(13),
      bin_end_notincl = paste(100),
      value = integrate(functions[[paste(pub_week + 3)]][[this_location]],
                         12.95, 15)$value
    ) %>%  
      mutate(value = ifelse(all(tail(wk3$value) == tail(wk3$value, n = 1)),
                            tail(wk3$value, n = 1), value)) %>%
      bind_rows(wk3, .)
    
    wk4 <- tibble(
      location = this_location,
      target = "4 wk ahead",
      type = "Bin",
      unit = "percent",
      bin_start_incl = paste(13),
      bin_end_notincl = paste(100),
      value = integrate(functions[[paste(pub_week + 4)]][[this_location]],
                         12.95, 15)$value
    ) %>%   
      mutate(value = ifelse(all(tail(wk4$value) == tail(wk4$value, n = 1)),
                            tail(wk4$value, n = 1), value)) %>%
      bind_rows(wk4, .)
    
    
    # Combine predictions together
    pred <- bind_rows(onset, pkwk, pkper, wk1, wk2, wk3, wk4, pred)
  }
  
  # Normalize probabilities, generate point forecasts, return prediction
  normalize_probs(pred) %>%
    bind_rows(generate_point_forecasts(.), .) %>%
    select(location, target, unit, type, bin_start_incl, 
           bin_end_notincl, value)

}
