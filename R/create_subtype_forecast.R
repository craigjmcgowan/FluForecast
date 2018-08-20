require(tidyverse)
require(FluSight)

# Create subtype weighted historical forecasts
create_subtype_forecast <- function(functions, virologic, pub_week, 
                                    season, prob_no_onset) {

  pred <- tibble()

  for(this_location in names(functions[["h1"]][["Season onset"]])) {
    
    temp_vir <- filter(virologic, location == this_location,
                       order_week == pub_week)
    
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
        value = (integrate(functions[["h1"]][["Season onset"]][[this_location]],
                           i - 0.5, i + 0.5)$value * temp_vir$h1per) +
          (integrate(functions[["h3"]][["Season onset"]][[this_location]],
                     i - 0.5, i + 0.5)$value * temp_vir$h3per)
      ) %>% bind_rows(onset, .)
      
      # Season peak week
      pkwk <- tibble(
        location = this_location,
        target = "Season peak week",
        type = "Bin",
        unit = "week",
        bin_start_incl = paste(i),
        bin_end_notincl = paste(i + 1),
        value = (integrate(functions[["h1"]][["Season peak week"]][[this_location]],
                           i - 0.5, i + 0.5)$value * temp_vir$h1per) +
          (integrate(functions[["h3"]][["Season peak week"]][[this_location]],
                     i - 0.5, i + 0.5)$value * temp_vir$h3per)
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
    prob_no <- prob_no_onset$h1_prob_no_onset[prob_no_onset$location == this_location] * temp_vir$h1per + 
      prob_no_onset$h3_prob_no_onset[prob_no_onset$location == this_location] * temp_vir$h3per
    
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
      value = (integrate(functions[["h1"]][["Season peak percentage"]][[this_location]],
                         0, 0.05)$value * temp_vir$h1per) +
        (integrate(functions[["h3"]][["Season peak percentage"]][[this_location]],
                   0, 0.05)$value * temp_vir$h3per)
    )
    
    wk1 <- tibble(
      location = this_location,
      target = "1 wk ahead",
      type = "Bin",
      unit = "percent",
      bin_start_incl = paste(0),
      bin_end_notincl = paste(0.1),
      value = (integrate(functions[["h1"]][[paste(pub_week + 1)]][[this_location]],
                         0, 0.05)$value * temp_vir$h1per) +
        (integrate(functions[["h3"]][[paste(pub_week + 1)]][[this_location]],
                   0, 0.05)$value * temp_vir$h3per)
    )
    
    wk2 <- tibble(
      location = this_location,
      target = "2 wk ahead",
      type = "Bin",
      unit = "percent",
      bin_start_incl = paste(0),
      bin_end_notincl = paste(0.1),
      value = (integrate(functions[["h1"]][[paste(pub_week + 2)]][[this_location]],
                         0, 0.05)$value * temp_vir$h1per) +
        (integrate(functions[["h3"]][[paste(pub_week + 2)]][[this_location]],
                   0, 0.05)$value * temp_vir$h3per)
    )
    
    wk3 <- tibble(
      location = this_location,
      target = "3 wk ahead",
      type = "Bin",
      unit = "percent",
      bin_start_incl = paste(0),
      bin_end_notincl = paste(0.1),
      value = (integrate(functions[["h1"]][[paste(pub_week + 3)]][[this_location]],
                         0, 0.05)$value * temp_vir$h1per) +
        (integrate(functions[["h3"]][[paste(pub_week + 3)]][[this_location]],
                   0, 0.05)$value * temp_vir$h3per)
    )
    
    wk4 <- tibble(
      location = this_location,
      target = "4 wk ahead",
      type = "Bin",
      unit = "percent",
      bin_start_incl = paste(0),
      bin_end_notincl = paste(0.1),
      value = (integrate(functions[["h1"]][[paste(pub_week + 4)]][[this_location]],
                         0, 0.05)$value * temp_vir$h1per) +
        (integrate(functions[["h3"]][[paste(pub_week + 4)]][[this_location]],
                   0, 0.05)$value * temp_vir$h3per)
    )
    
    
    for (i in seq(0.1, 12.9, 0.1)) {
      
      pkper <- tibble(
        location = this_location,
        target = "Season peak percentage",
        type = "Bin",
        unit = "percent",
        bin_start_incl = paste(i),
        bin_end_notincl = paste(i + 0.1),
        value = (integrate(functions[["h1"]][["Season peak percentage"]][[this_location]],
                           i - 0.05, i + 0.05)$value * temp_vir$h1per) +
          (integrate(functions[["h3"]][["Season peak percentage"]][[this_location]],
                     i - 0.05, i + 0.05)$value * temp_vir$h3per)
      ) %>% bind_rows(pkper, .)
      
      wk1 <- tibble(
        location = this_location,
        target = "1 wk ahead",
        type = "Bin",
        unit = "percent",
        bin_start_incl = paste(i),
        bin_end_notincl = paste(i + 0.1),
        value = (integrate(functions[["h1"]][[paste(pub_week + 1)]][[this_location]],
                           i - 0.05, i + 0.05)$value * temp_vir$h1per) +
          (integrate(functions[["h3"]][[paste(pub_week + 1)]][[this_location]],
                     i - 0.05, i + 0.05)$value * temp_vir$h3per)
      ) %>% bind_rows(wk1, .)
      
      wk2 <- tibble(
        location = this_location,
        target = "2 wk ahead",
        type = "Bin",
        unit = "percent",
        bin_start_incl = paste(i),
        bin_end_notincl = paste(i + 0.1),
        value = (integrate(functions[["h1"]][[paste(pub_week + 2)]][[this_location]],
                           i - 0.05, i + 0.05)$value * temp_vir$h1per) +
          (integrate(functions[["h3"]][[paste(pub_week + 2)]][[this_location]],
                     i - 0.05, i + 0.05)$value * temp_vir$h3per)
      ) %>% bind_rows(wk2, .)
      
      wk3 <- tibble(
        location = this_location,
        target = "3 wk ahead",
        type = "Bin",
        unit = "percent",
        bin_start_incl = paste(i),
        bin_end_notincl = paste(i + 0.1),
        value = (integrate(functions[["h1"]][[paste(pub_week + 3)]][[this_location]],
                           i - 0.05, i + 0.05)$value * temp_vir$h1per) +
          (integrate(functions[["h3"]][[paste(pub_week + 3)]][[this_location]],
                     i - 0.05, i + 0.05)$value * temp_vir$h3per)
      ) %>% bind_rows(wk3, .)
      
      wk4 <- tibble(
        location = this_location,
        target = "4 wk ahead",
        type = "Bin",
        unit = "percent",
        bin_start_incl = paste(i),
        bin_end_notincl = paste(i + 0.1),
        value = (integrate(functions[["h1"]][[paste(pub_week + 4)]][[this_location]],
                           i - 0.05, i + 0.05)$value * temp_vir$h1per) +
          (integrate(functions[["h3"]][[paste(pub_week + 4)]][[this_location]],
                     i - 0.05, i + 0.05)$value * temp_vir$h3per)
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
      value = (integrate(functions[["h1"]][["Season peak percentage"]][[this_location]],
                         12.95, 15)$value * temp_vir$h1per) +
        (integrate(functions[["h3"]][["Season peak percentage"]][[this_location]],
                   12.95, 15)$value * temp_vir$h3per)
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
      value = (integrate(functions[["h1"]][[paste(pub_week + 1)]][[this_location]],
                         12.95, 15)$value * temp_vir$h1per) +
        (integrate(functions[["h3"]][[paste(pub_week + 1)]][[this_location]],
                   12.95, 15)$value * temp_vir$h3per)
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
      value = (integrate(functions[["h1"]][[paste(pub_week + 2)]][[this_location]],
                         12.95, 15)$value * temp_vir$h1per) +
        (integrate(functions[["h3"]][[paste(pub_week + 2)]][[this_location]],
                   12.95, 15)$value * temp_vir$h3per)
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
      value = (integrate(functions[["h1"]][[paste(pub_week + 3)]][[this_location]],
                         12.95, 15)$value * temp_vir$h1per) +
        (integrate(functions[["h3"]][[paste(pub_week + 3)]][[this_location]],
                   12.95, 15)$value * temp_vir$h3per)
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
      value = (integrate(functions[["h1"]][[paste(pub_week + 4)]][[this_location]],
                         12.95, 15)$value * temp_vir$h1per) +
        (integrate(functions[["h3"]][[paste(pub_week + 4)]][[this_location]],
                   12.95, 15)$value * temp_vir$h3per)
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
