# Create subtype-specific densities
# Make densities based on provided data -----
create_subtype_densities <- function(ili_df, vir_ssn_per) {
  
  # Only keep data from 2000/2001 season on, except 2009/2010
  train_ili <- filter(ili_df, season != "2009/2010",
                      year >= 2000, season != "1999/2000") %>%
    select(location, season, week, year, ILI) %>%
    mutate(order_week = week_inorder(week, season))
  
  # Generate past onsets from years with baseline data
  train_onsets <- train_ili %>%
    filter(year >= 2007, season != "2006/2007",
           (week >= 40 | week <= 20)) %>%
    group_by(location, season) %>%
    do(create_onset(., region = .$location[1], 
                    year = as.numeric(substr(.$season[1], 1, 4)))) %>%
    ungroup()
  
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
        TRUE ~ as.numeric(bin_start_incl))) %>%
      left_join(vir_ssn_per, by = c("season", "location"))
    
    densities[["h1"]][["Season onset"]][[this_location]] <- 
      density(temp_onset$bin_start_incl, bw = "SJ",
              weights = temp_onset$h1per / sum(temp_onset$h1per))
    densities[["h3"]][["Season onset"]][[this_location]] <- 
      density(temp_onset$bin_start_incl, bw = "SJ",
              weights = temp_onset$h3per / sum(temp_onset$h3per))
    
    # Peak week densities
    temp_peak_wk <- train_peak %>%
      filter(location == this_location, target == "Season peak week") %>%
      mutate(bin_start_incl = case_when(
        as.numeric(bin_start_incl) < 40 & 
          season %in% c("1997/1998", "2003/2004", "2008/2009", "2014/2015") ~
          as.numeric(bin_start_incl) + 53,
        as.numeric(bin_start_incl) < 40 ~ 
          as.numeric(bin_start_incl) + 52,
        TRUE ~ as.numeric(bin_start_incl))) %>%
      left_join(vir_ssn_per, by = c("season", "location"))
    
    densities[["h1"]][["Season peak week"]][[this_location]] <- 
      density(temp_peak_wk$bin_start_incl, bw = "SJ",
              weights = temp_peak_wk$h1per / sum(temp_peak_wk$h1per))
    densities[["h3"]][["Season peak week"]][[this_location]] <- 
      density(temp_peak_wk$bin_start_incl, bw = "SJ",
              weights = temp_peak_wk$h3per / sum(temp_peak_wk$h3per))
    
    # Peak percentage densities
    temp_peak_per <- train_peak %>%
      mutate(bin_start_incl = as.numeric(bin_start_incl)) %>%
      filter(location == this_location, target == "Season peak percentage") %>%
      left_join(vir_ssn_per, by = c("season", "location"))
    
    densities[["h1"]][["Season peak percentage"]][[this_location]] <- 
      density(temp_peak_per$bin_start_incl, bw = "SJ",
              weights = temp_peak_per$h1per / sum(temp_peak_per$h1per))
    densities[["h3"]][["Season peak percentage"]][[this_location]] <- 
      density(temp_peak_per$bin_start_incl, bw = "SJ",
              weights = temp_peak_per$h3per / sum(temp_peak_per$h3per))
    
    
    # Densities for each week of the season
    for(this_week in 40:77) {
      temp_ili <- filter(train_ili, location == this_location,
                         order_week == this_week) %>%
        left_join(vir_ssn_per, by = c("season", "location"))
      
      densities[["h1"]][[paste(this_week)]][[this_location]] <-
        density(temp_ili$ILI, bw = "SJ",
                weights = temp_ili$h1per / sum(temp_ili$h1per))
      densities[["h3"]][[paste(this_week)]][[this_location]] <-
        density(temp_ili$ILI, bw = "SJ",
                weights = temp_ili$h3per / sum(temp_ili$h3per))
    }
  }
  
  return(densities)
  
}


