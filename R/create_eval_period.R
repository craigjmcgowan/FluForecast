# Create evaluation period for particular year

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
