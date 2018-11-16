# Pull in and save datasets from EpiData and other sources 
library(tidyverse)
library(cdcfluview)
library(MMWRweek)
library(gtrendsR)
library(lubridate)
library(FluSight)
# library(epiforecast)

# Load functions
source("R/utils.R")
source('R/EpiDataAPI.R')

# Save current MMWR week in format Epidata wants
current_week <- as.numeric(paste0(MMWRweek(Sys.Date())[[1]], MMWRweek(Sys.Date())[[2]] ))

pull_week <- case_when(
  MMWRweek(Sys.Date())[3] == 7 & substr(current_week, 5, 6) == "01" ~ 
    as.numeric(paste0(as.numeric(substr(current_week, 1, 4)) - 1, "52")),
  MMWRweek(Sys.Date())[3] == 7 ~
    current_week - 1,
  substr(current_week, 5, 6) == "01" ~ 
    as.numeric(paste0(as.numeric(substr(current_week, 1, 4)) - 1, "51")),
  substr(current_week, 5, 6) == "02" ~ 
    as.numeric(paste0(as.numeric(substr(current_week, 1, 4)) - 1, "52")),
  TRUE ~ current_week - 2
)

# Fetch wILI data from EpiData API -------
# Save ILINet values for previous 26 weeks for each week's publication available
ili_init_pub_list <- list()

for(i in c(201040:201052, 201101:201152, 201201:201252, 201301:201338,
           201340:201352, 201401:201453, 201501:201552, 201601:201652,
           201701:201752, 201801:pull_week)) {
  ili_init_pub_list[[paste(i)]] <- pull_initpub_epidata(i) %>%
    mutate(year = as.integer(substr(epiweek, 1, 4)),
           week = as.integer(substr(epiweek, 5, 6)),
           season = ifelse(week >= 40,
                           paste0(year, "/", year + 1),
                           paste0(year - 1, "/", year)),
           location = case_when(
             region == "nat" ~ "US National",
             region == "hhs1" ~ "HHS Region 1",
             region == "hhs2" ~ "HHS Region 2",
             region == "hhs3" ~ "HHS Region 3",
             region == "hhs4" ~ "HHS Region 4",
             region == "hhs5" ~ "HHS Region 5",
             region == "hhs6" ~ "HHS Region 6",
             region == "hhs7" ~ "HHS Region 7",
             region == "hhs8" ~ "HHS Region 8",
             region == "hhs9" ~ "HHS Region 9",
             region == "hhs10" ~ "HHS Region 10"
           )) %>%
    rename(ILI = wili) %>%
    select(-ili) %>%
    # Make ILI always > 0
    mutate(ILI = case_when(near(ILI, 0) ~ 0.1,
                           TRUE ~ ILI))
}

# Original ILI data for each week
ili_orig <- bind_rows(ili_init_pub_list) %>%
  filter(lag == 0)
  
# ILI currently ------
ili_current <- bind_rows(
  pull_curr_epidata(199740, 200139),
  pull_curr_epidata(200140, 200439),
  pull_curr_epidata(200440, 200739),
  pull_curr_epidata(200740, 201039),
  pull_curr_epidata(201040, 201339),
  pull_curr_epidata(201340, 201639),
  pull_curr_epidata(201640, pull_week)) %>%
  mutate(year = as.integer(substr(epiweek, 1, 4)),
         week = as.integer(substr(epiweek, 5, 6)),
         season = ifelse(week >= 40,
                         paste0(year, "/", year + 1),
                         paste0(year - 1, "/", year)),
         location = case_when(
           region == "nat" ~ "US National",
           region == "hhs1" ~ "HHS Region 1",
             region == "hhs2" ~ "HHS Region 2",
             region == "hhs3" ~ "HHS Region 3",
             region == "hhs4" ~ "HHS Region 4",
             region == "hhs5" ~ "HHS Region 5",
             region == "hhs6" ~ "HHS Region 6",
             region == "hhs7" ~ "HHS Region 7",
             region == "hhs8" ~ "HHS Region 8",
             region == "hhs9" ~ "HHS Region 9",
             region == "hhs10" ~ "HHS Region 10"
         )) %>%
  rename(ILI = wili) %>%
  select(-ili) %>%
  # Make ILI always > 0
  mutate(ILI = case_when(near(ILI, 0) ~ 0.1,
                         TRUE ~ ILI))

# Create measure of backfill
ili_backfill <- left_join(
  select(ili_orig, location, season, year, week, orig_ILI = ILI),
  select(ili_current, location, season, year, week, final_ILI = ILI),
  by = c("location", "season", "year", "week")
) %>%
  mutate(backfill = final_ILI - orig_ILI)

ili_backfill_avg <- ili_backfill %>%
  group_by(location, week) %>%
  summarize(avg_backfill = mean(backfill),
            sd_backfill = sd(backfill))

save(ili_orig, ili_current, ili_init_pub_list, ili_backfill, ili_backfill_avg,
     file = "Data/ili.Rdata")

# Fetch virologic data from CDC -----
virologic_national <- who_nrevss(region = "national", years = c(1997:2018))
virologic_region <- who_nrevss(region = "hhs", years = c(1997:2018))

virologic_before_1516 <- bind_rows(
  virologic_national[[1]], virologic_region[[1]]
) %>%
  mutate_at(c("a_2009_h1n1", "a_h1", "a_h3", "a_subtyping_not_performed",
              "a_unable_to_subtype", "b", "h3n2v"),
            as.integer)

virologic_ph_lab <- bind_rows(
  virologic_national[[2]], virologic_region[[2]]
)

virologic_combined <- bind_rows(
  virologic_before_1516, virologic_ph_lab
) %>%
  # Create season indicator
  mutate(season = ifelse(week >= 40,
                         paste0(year, "/", year + 1),
                         paste0(year - 1, "/", year)),
         region = ifelse(region == "National", "US National", 
                         paste("HHS", region))) %>%
  # Create mock H1 and H3 percents if all A samples had been tested
  rowwise() %>%
  mutate(pos_samples = sum(a_h1, a_2009_h1n1, a_h3, a_subtyping_not_performed, 
           a_unable_to_subtype, b, bvic, byam, na.rm = TRUE),
         a_pos_samples = sum(a_h1, a_2009_h1n1, a_h3, a_subtyping_not_performed,
                             a_unable_to_subtype, na.rm = TRUE),
         b_pos_samples = sum(b, bvic, byam, na.rm = TRUE)) %>% #,
         # h1_per_samples = (sum(a_h1, a_2009_h1n1, na.rm = TRUE) / 
         #                     sum(a_h1, a_2009_h1n1, a_h3, na.rm = TRUE)) *
         #   sum(a_h1, a_2009_h1n1, a_h3, a_subtyping_not_performed, 
         #       a_unable_to_subtype, na.rm = T) / pos_samples,
         # h3_per_samples = a_h3 / sum(a_h1, a_2009_h1n1, a_h3, na.rm = TRUE) *
         #   sum(a_h1, a_2009_h1n1, a_h3, a_subtyping_not_performed, 
         #       a_unable_to_subtype, na.rm = T) / pos_samples,
         # b_per_samples = b_pos_samples / pos_samples) %>%
  ungroup() %>%
  mutate_at(vars(c("a_2009_h1n1", "a_h1", "a_h3", "a_subtyping_not_performed",
                   # "h1_per_samples", "h3_per_samples", "b_per_samples",
                   "a_unable_to_subtype", "b", "bvic", "byam")),
            function(x) ifelse(is.na(x), 0, x)) %>%
  # Create cumulative influenza percentage measures
  group_by(season, region) %>%
  arrange(wk_date, .by_group = T) %>%
  mutate(h1sum = cumsum(a_2009_h1n1) + cumsum(a_h1),
         h3sum = cumsum(a_h3),
         asum = cumsum(a_pos_samples),
         bsum = cumsum(b_pos_samples),
         # Subtype percentages of typed As
         h1per_of_a = ifelse(is.na(h1sum / (h1sum + h3sum)), 0.5, 
                             h1sum / (h1sum + h3sum)),
         h3per_of_a = ifelse(is.na(h3sum / (h1sum + h3sum)), 0.5, 
                             h3sum / (h1sum + h3sum)),
         # Cumulative percentages of each type, assuming A type %s are representative
         cum_bper = bsum / (asum + bsum),
         cum_h1per = h1per_of_a * asum / (asum + bsum),
         cum_h3per = h3per_of_a * asum / (asum + bsum)) %>%
  # If no samples reported, make each type 33%
  mutate(cum_h1per = ifelse(is.na(cum_h1per), 1/3, cum_h1per),
         cum_h3per = ifelse(is.na(cum_h3per), 1/3, cum_h3per),
         cum_bper = ifelse(is.na(cum_bper), 1/3, cum_bper)) %>% 
  ungroup() %>%
  select(location = region, season, year, week, a_h1, a_2009_h1n1, a_h3, 
         a_subtyping_not_performed, a_unable_to_subtype, b, bvic, byam, 
         h1per_of_a, h3per_of_a, cum_h1per, cum_h3per, cum_bper)

save(virologic_combined,
     file = "Data/virologic.Rdata")


