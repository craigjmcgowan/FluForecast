# Google Trends data
library(tidyverse)
library(MMWRweek)
library(gtrendsR)
library(lubridate)

load("Data/Gtrends.Rdata")
load("Data/virologic.Rdata")
load("Data/ili.Rdata")
source("R/utils.R")

# Explore trends by state - 1 state in each region ------------

# Region 1 - Massachusetts
MA_flu_0407 <- gtrends(keyword = "flu",
                       geo = "US-MA",
                       time = paste(MMWRweek2Date(2004, 41), MMWRweek2Date(2007, 40)))$interest_over_time %>%
  mutate(hits = case_when(hits == "<1" ~ 1,
                          TRUE ~ as.numeric(hits)))

MA_flu_0611 <- gtrends(keyword = "flu",
                       geo = "US-MA",
                       time = paste(MMWRweek2Date(2006, 41), MMWRweek2Date(2011, 40)))$interest_over_time %>%
  mutate(hits = case_when(hits == "<1" ~ 1,
                          TRUE ~ as.numeric(hits)))

MA_flu_1015 <- gtrends(keyword = "flu",
                       geo = "US-MA",
                       time = paste(MMWRweek2Date(2010, 41), MMWRweek2Date(2015, 40)))$interest_over_time %>%
  mutate(hits = case_when(hits == "<1" ~ 1,
                          TRUE ~ as.numeric(hits)))

MA_flu_1419 <- gtrends(keyword = "flu",
                       geo = "US-MA")$interest_over_time %>%
  mutate(hits = case_when(hits == "<1" ~ 1,
                          TRUE ~ as.numeric(hits)))

# Set up ratios to normalize all Gtrends data to scale from 2010-2015
MA_gratio_0607 <- US_flu_ratio(MA_flu_0407, MA_flu_0611)
MA_gratio_1011 <- US_flu_ratio(MA_flu_0611, MA_flu_1015)
MA_inv_gratio_1415 <- US_flu_ratio(MA_flu_1419, MA_flu_1015)

# Merge Gtrends data and rescale to 2010-2015 scale
gtrend_MA_flu_merge <- filter(MA_flu_0407, !date %in% MA_flu_0611$date) %>%
  mutate(hits = hits / MA_gratio_0607) %>%
  bind_rows(MA_flu_0611) %>%
  filter(!date %in% MA_flu_1015$date) %>%
  mutate(hits = hits / MA_gratio_1011) %>%
  bind_rows(MA_flu_1015) %>%
  filter(!date %in% MA_flu_1419$date) %>%
  bind_rows(MA_flu_1419 %>%
              mutate(hits = hits / MA_inv_gratio_1415)) %>%
  select(date, hits) %>%
  mutate(week = MMWRweek(date)[[2]],
         year = MMWRweek(date)[[1]],
         season = ifelse(week >= 40,
                         paste0(year, "/", year + 1),
                         paste0(year - 1, "/", year)),
         hits = case_when(near(hits, 0) ~ 1,
                          TRUE ~ hits))

# Region 2 - New York
NY_flu_0407 <- gtrends(keyword = "flu",
                       geo = "US-NY",
                       time = paste(MMWRweek2Date(2004, 41), MMWRweek2Date(2007, 40)))$interest_over_time %>%
  mutate(hits = case_when(hits == "<1" ~ 1,
                          TRUE ~ as.numeric(hits)))

NY_flu_0611 <- gtrends(keyword = "flu",
                       geo = "US-NY",
                       time = paste(MMWRweek2Date(2006, 41), MMWRweek2Date(2011, 40)))$interest_over_time %>%
  mutate(hits = case_when(hits == "<1" ~ 1,
                          TRUE ~ as.numeric(hits)))

NY_flu_1015 <- gtrends(keyword = "flu",
                       geo = "US-NY",
                       time = paste(MMWRweek2Date(2010, 41), MMWRweek2Date(2015, 40)))$interest_over_time %>%
  mutate(hits = case_when(hits == "<1" ~ 1,
                          TRUE ~ as.numeric(hits)))

NY_flu_1419 <- gtrends(keyword = "flu",
                       geo = "US-NY")$interest_over_time %>%
  mutate(hits = case_when(hits == "<1" ~ 1,
                          TRUE ~ as.numeric(hits)))

# Set up ratios to normalize all Gtrends data to scale from 2010-2015
NY_gratio_0607 <- US_flu_ratio(NY_flu_0407, NY_flu_0611)
NY_gratio_1011 <- US_flu_ratio(NY_flu_0611, NY_flu_1015)
NY_inv_gratio_1415 <- US_flu_ratio(NY_flu_1419, NY_flu_1015)

# Merge Gtrends data and rescale to 2010-2015 scale
gtrend_NY_flu_merge <- filter(NY_flu_0407, !date %in% NY_flu_0611$date) %>%
  mutate(hits = hits / NY_gratio_0607) %>%
  bind_rows(NY_flu_0611) %>%
  filter(!date %in% NY_flu_1015$date) %>%
  mutate(hits = hits / NY_gratio_1011) %>%
  bind_rows(NY_flu_1015) %>%
  filter(!date %in% NY_flu_1419$date) %>%
  bind_rows(NY_flu_1419 %>%
              mutate(hits = hits / NY_inv_gratio_1415)) %>%
  select(date, hits) %>%
  mutate(week = MMWRweek(date)[[2]],
         year = MMWRweek(date)[[1]],
         season = ifelse(week >= 40,
                         paste0(year, "/", year + 1),
                         paste0(year - 1, "/", year)),
         hits = case_when(near(hits, 0) ~ 1,
                          TRUE ~ hits))

# Region 3 - Pennsylvania
PA_flu_0407 <- gtrends(keyword = "flu",
                       geo = "US-PA",
                       time = paste(MMWRweek2Date(2004, 41), MMWRweek2Date(2007, 40)))$interest_over_time %>%
  mutate(hits = case_when(hits == "<1" ~ 1,
                          TRUE ~ as.numeric(hits)))

PA_flu_0611 <- gtrends(keyword = "flu",
                       geo = "US-PA",
                       time = paste(MMWRweek2Date(2006, 41), MMWRweek2Date(2011, 40)))$interest_over_time %>%
  mutate(hits = case_when(hits == "<1" ~ 1,
                          TRUE ~ as.numeric(hits)))

PA_flu_1015 <- gtrends(keyword = "flu",
                       geo = "US-PA",
                       time = paste(MMWRweek2Date(2010, 41), MMWRweek2Date(2015, 40)))$interest_over_time %>%
  mutate(hits = case_when(hits == "<1" ~ 1,
                          TRUE ~ as.numeric(hits)))

PA_flu_1419 <- gtrends(keyword = "flu",
                       geo = "US-PA")$interest_over_time %>%
  mutate(hits = case_when(hits == "<1" ~ 1,
                          TRUE ~ as.numeric(hits)))

# Set up ratios to normalize all Gtrends data to scale from 2010-2015
PA_gratio_0607 <- US_flu_ratio(PA_flu_0407, PA_flu_0611)
PA_gratio_1011 <- US_flu_ratio(PA_flu_0611, PA_flu_1015)
PA_inv_gratio_1415 <- US_flu_ratio(PA_flu_1419, PA_flu_1015)

# Merge Gtrends data and rescale to 2010-2015 scale
gtrend_PA_flu_merge <- filter(PA_flu_0407, !date %in% PA_flu_0611$date) %>%
  mutate(hits = hits / PA_gratio_0607) %>%
  bind_rows(PA_flu_0611) %>%
  filter(!date %in% PA_flu_1015$date) %>%
  mutate(hits = hits / PA_gratio_1011) %>%
  bind_rows(PA_flu_1015) %>%
  filter(!date %in% PA_flu_1419$date) %>%
  bind_rows(PA_flu_1419 %>%
              mutate(hits = hits / PA_inv_gratio_1415)) %>%
  select(date, hits) %>%
  mutate(week = MMWRweek(date)[[2]],
         year = MMWRweek(date)[[1]],
         season = ifelse(week >= 40,
                         paste0(year, "/", year + 1),
                         paste0(year - 1, "/", year)),
         hits = case_when(near(hits, 0) ~ 1,
                          TRUE ~ hits))


# Region 4 - Florida
FL_flu_0407 <- gtrends(keyword = "flu",
                       geo = "US-FL",
                       time = paste(MMWRweek2Date(2004, 41), MMWRweek2Date(2007, 40)))$interest_over_time %>%
  mutate(hits = case_when(hits == "<1" ~ 1,
                          TRUE ~ as.numeric(hits)))

FL_flu_0611 <- gtrends(keyword = "flu",
                       geo = "US-FL",
                       time = paste(MMWRweek2Date(2006, 41), MMWRweek2Date(2011, 40)))$interest_over_time %>%
  mutate(hits = case_when(hits == "<1" ~ 1,
                          TRUE ~ as.numeric(hits)))

FL_flu_1015 <- gtrends(keyword = "flu",
                       geo = "US-FL",
                       time = paste(MMWRweek2Date(2010, 41), MMWRweek2Date(2015, 40)))$interest_over_time %>%
  mutate(hits = case_when(hits == "<1" ~ 1,
                          TRUE ~ as.numeric(hits)))

FL_flu_1419 <- gtrends(keyword = "flu",
                       geo = "US-FL")$interest_over_time %>%
  mutate(hits = case_when(hits == "<1" ~ 1,
                          TRUE ~ as.numeric(hits)))

# Set up ratios to normalize all Gtrends data to scale from 2010-2015
FL_gratio_0607 <- US_flu_ratio(FL_flu_0407, FL_flu_0611)
FL_gratio_1011 <- US_flu_ratio(FL_flu_0611, FL_flu_1015)
FL_inv_gratio_1415 <- US_flu_ratio(FL_flu_1419, FL_flu_1015)

# Merge Gtrends data and rescale to 2010-2015 scale
gtrend_FL_flu_merge <- filter(FL_flu_0407, !date %in% FL_flu_0611$date) %>%
  mutate(hits = hits / FL_gratio_0607) %>%
  bind_rows(FL_flu_0611) %>%
  filter(!date %in% FL_flu_1015$date) %>%
  mutate(hits = hits / FL_gratio_1011) %>%
  bind_rows(FL_flu_1015) %>%
  filter(!date %in% FL_flu_1419$date) %>%
  bind_rows(FL_flu_1419 %>%
              mutate(hits = hits / FL_inv_gratio_1415)) %>%
  select(date, hits) %>%
  mutate(week = MMWRweek(date)[[2]],
         year = MMWRweek(date)[[1]],
         season = ifelse(week >= 40,
                         paste0(year, "/", year + 1),
                         paste0(year - 1, "/", year)),
         hits = case_when(near(hits, 0) ~ 1,
                          TRUE ~ hits))

# Region 5 - Illinois
IL_flu_0407 <- gtrends(keyword = "flu",
                       geo = "US-IL",
                       time = paste(MMWRweek2Date(2004, 41), MMWRweek2Date(2007, 40)))$interest_over_time %>%
  mutate(hits = case_when(hits == "<1" ~ 1,
                          TRUE ~ as.numeric(hits)))

IL_flu_0611 <- gtrends(keyword = "flu",
                       geo = "US-IL",
                       time = paste(MMWRweek2Date(2006, 41), MMWRweek2Date(2011, 40)))$interest_over_time %>%
  mutate(hits = case_when(hits == "<1" ~ 1,
                          TRUE ~ as.numeric(hits)))

IL_flu_1015 <- gtrends(keyword = "flu",
                       geo = "US-IL",
                       time = paste(MMWRweek2Date(2010, 41), MMWRweek2Date(2015, 40)))$interest_over_time %>%
  mutate(hits = case_when(hits == "<1" ~ 1,
                          TRUE ~ as.numeric(hits)))

IL_flu_1419 <- gtrends(keyword = "flu",
                       geo = "US-IL")$interest_over_time %>%
  mutate(hits = case_when(hits == "<1" ~ 1,
                          TRUE ~ as.numeric(hits)))

# Set up ratios to normalize all Gtrends data to scale from 2010-2015
IL_gratio_0607 <- US_flu_ratio(IL_flu_0407, IL_flu_0611)
IL_gratio_1011 <- US_flu_ratio(IL_flu_0611, IL_flu_1015)
IL_inv_gratio_1415 <- US_flu_ratio(IL_flu_1419, IL_flu_1015)

# Merge Gtrends data and rescale to 2010-2015 scale
gtrend_IL_flu_merge <- filter(IL_flu_0407, !date %in% IL_flu_0611$date) %>%
  mutate(hits = hits / IL_gratio_0607) %>%
  bind_rows(IL_flu_0611) %>%
  filter(!date %in% IL_flu_1015$date) %>%
  mutate(hits = hits / IL_gratio_1011) %>%
  bind_rows(IL_flu_1015) %>%
  filter(!date %in% IL_flu_1419$date) %>%
  bind_rows(IL_flu_1419 %>%
              mutate(hits = hits / IL_inv_gratio_1415)) %>%
  select(date, hits) %>%
  mutate(week = MMWRweek(date)[[2]],
         year = MMWRweek(date)[[1]],
         season = ifelse(week >= 40,
                         paste0(year, "/", year + 1),
                         paste0(year - 1, "/", year)),
         hits = case_when(near(hits, 0) ~ 1,
                          TRUE ~ hits))

# Region 6 - Texas
TX_flu_0407 <- gtrends(keyword = "flu",
                       geo = "US-TX",
                       time = paste(MMWRweek2Date(2004, 41), MMWRweek2Date(2007, 40)))$interest_over_time %>%
  mutate(hits = case_when(hits == "<1" ~ 1,
                          TRUE ~ as.numeric(hits)))

TX_flu_0611 <- gtrends(keyword = "flu",
                       geo = "US-TX",
                       time = paste(MMWRweek2Date(2006, 41), MMWRweek2Date(2011, 40)))$interest_over_time %>%
  mutate(hits = case_when(hits == "<1" ~ 1,
                          TRUE ~ as.numeric(hits)))

TX_flu_1015 <- gtrends(keyword = "flu",
                       geo = "US-TX",
                       time = paste(MMWRweek2Date(2010, 41), MMWRweek2Date(2015, 40)))$interest_over_time %>%
  mutate(hits = case_when(hits == "<1" ~ 1,
                          TRUE ~ as.numeric(hits)))

TX_flu_1419 <- gtrends(keyword = "flu",
                       geo = "US-TX")$interest_over_time %>%
  mutate(hits = case_when(hits == "<1" ~ 1,
                          TRUE ~ as.numeric(hits)))

# Set up ratios to normalize all Gtrends data to scale from 2010-2015
TX_gratio_0607 <- US_flu_ratio(TX_flu_0407, TX_flu_0611)
TX_gratio_1011 <- US_flu_ratio(TX_flu_0611, TX_flu_1015)
TX_inv_gratio_1415 <- US_flu_ratio(TX_flu_1419, TX_flu_1015)

# Merge Gtrends data and rescale to 2010-2015 scale
gtrend_TX_flu_merge <- filter(TX_flu_0407, !date %in% TX_flu_0611$date) %>%
  mutate(hits = hits / TX_gratio_0607) %>%
  bind_rows(TX_flu_0611) %>%
  filter(!date %in% TX_flu_1015$date) %>%
  mutate(hits = hits / TX_gratio_1011) %>%
  bind_rows(TX_flu_1015) %>%
  filter(!date %in% TX_flu_1419$date) %>%
  bind_rows(TX_flu_1419 %>%
              mutate(hits = hits / TX_inv_gratio_1415)) %>%
  select(date, hits) %>%
  mutate(week = MMWRweek(date)[[2]],
         year = MMWRweek(date)[[1]],
         season = ifelse(week >= 40,
                         paste0(year, "/", year + 1),
                         paste0(year - 1, "/", year)),
         hits = case_when(near(hits, 0) ~ 1,
                          TRUE ~ hits))

# Region 7 - Missouri
MO_flu_0407 <- gtrends(keyword = "flu",
                       geo = "US-MO",
                       time = paste(MMWRweek2Date(2004, 41), MMWRweek2Date(2007, 40)))$interest_over_time %>%
  mutate(hits = case_when(hits == "<1" ~ 1,
                          TRUE ~ as.numeric(hits)))

MO_flu_0611 <- gtrends(keyword = "flu",
                       geo = "US-MO",
                       time = paste(MMWRweek2Date(2006, 41), MMWRweek2Date(2011, 40)))$interest_over_time %>%
  mutate(hits = case_when(hits == "<1" ~ 1,
                          TRUE ~ as.numeric(hits)))

MO_flu_1015 <- gtrends(keyword = "flu",
                       geo = "US-MO",
                       time = paste(MMWRweek2Date(2010, 41), MMWRweek2Date(2015, 40)))$interest_over_time %>%
  mutate(hits = case_when(hits == "<1" ~ 1,
                          TRUE ~ as.numeric(hits)))

MO_flu_1419 <- gtrends(keyword = "flu",
                       geo = "US-MO")$interest_over_time %>%
  mutate(hits = case_when(hits == "<1" ~ 1,
                          TRUE ~ as.numeric(hits)))

# Set up ratios to normalize all Gtrends data to scale from 2010-2015
MO_gratio_0607 <- US_flu_ratio(MO_flu_0407, MO_flu_0611)
MO_gratio_1011 <- US_flu_ratio(MO_flu_0611, MO_flu_1015)
MO_inv_gratio_1415 <- US_flu_ratio(MO_flu_1419, MO_flu_1015)

# Merge Gtrends data and rescale to 2010-2015 scale
gtrend_MO_flu_merge <- filter(MO_flu_0407, !date %in% MO_flu_0611$date) %>%
  mutate(hits = hits / MO_gratio_0607) %>%
  bind_rows(MO_flu_0611) %>%
  filter(!date %in% MO_flu_1015$date) %>%
  mutate(hits = hits / MO_gratio_1011) %>%
  bind_rows(MO_flu_1015) %>%
  filter(!date %in% MO_flu_1419$date) %>%
  bind_rows(MO_flu_1419 %>%
              mutate(hits = hits / MO_inv_gratio_1415)) %>%
  select(date, hits) %>%
  mutate(week = MMWRweek(date)[[2]],
         year = MMWRweek(date)[[1]],
         season = ifelse(week >= 40,
                         paste0(year, "/", year + 1),
                         paste0(year - 1, "/", year)),
         hits = case_when(near(hits, 0) ~ 1,
                          TRUE ~ hits))

# Region 8 - Colorado
CO_flu_0407 <- gtrends(keyword = "flu",
                       geo = "US-CO",
                       time = paste(MMWRweek2Date(2004, 41), MMWRweek2Date(2007, 40)))$interest_over_time %>%
  mutate(hits = case_when(hits == "<1" ~ 1,
                          TRUE ~ as.numeric(hits)))

CO_flu_0611 <- gtrends(keyword = "flu",
                       geo = "US-CO",
                       time = paste(MMWRweek2Date(2006, 41), MMWRweek2Date(2011, 40)))$interest_over_time %>%
  mutate(hits = case_when(hits == "<1" ~ 1,
                          TRUE ~ as.numeric(hits)))

CO_flu_1015 <- gtrends(keyword = "flu",
                       geo = "US-CO",
                       time = paste(MMWRweek2Date(2010, 41), MMWRweek2Date(2015, 40)))$interest_over_time %>%
  mutate(hits = case_when(hits == "<1" ~ 1,
                          TRUE ~ as.numeric(hits)))

CO_flu_1419 <- gtrends(keyword = "flu",
                       geo = "US-CO")$interest_over_time %>%
  mutate(hits = case_when(hits == "<1" ~ 1,
                          TRUE ~ as.numeric(hits)))

# Set up ratios to normalize all Gtrends data to scale from 2010-2015
CO_gratio_0607 <- US_flu_ratio(CO_flu_0407, CO_flu_0611)
CO_gratio_1011 <- US_flu_ratio(CO_flu_0611, CO_flu_1015)
CO_inv_gratio_1415 <- US_flu_ratio(CO_flu_1419, CO_flu_1015)

# Merge Gtrends data and rescale to 2010-2015 scale
gtrend_CO_flu_merge <- filter(CO_flu_0407, !date %in% CO_flu_0611$date) %>%
  mutate(hits = hits / CO_gratio_0607) %>%
  bind_rows(CO_flu_0611) %>%
  filter(!date %in% CO_flu_1015$date) %>%
  mutate(hits = hits / CO_gratio_1011) %>%
  bind_rows(CO_flu_1015) %>%
  filter(!date %in% CO_flu_1419$date) %>%
  bind_rows(CO_flu_1419 %>%
              mutate(hits = hits / CO_inv_gratio_1415)) %>%
  select(date, hits) %>%
  mutate(week = MMWRweek(date)[[2]],
         year = MMWRweek(date)[[1]],
         season = ifelse(week >= 40,
                         paste0(year, "/", year + 1),
                         paste0(year - 1, "/", year)),
         hits = case_when(near(hits, 0) ~ 1,
                          TRUE ~ hits))

# Region 9 - California
CA_flu_0407 <- gtrends(keyword = "flu",
                       geo = "US-CA",
                       time = paste(MMWRweek2Date(2004, 41), MMWRweek2Date(2007, 40)))$interest_over_time %>%
  mutate(hits = case_when(hits == "<1" ~ 1,
                          TRUE ~ as.numeric(hits)))

CA_flu_0611 <- gtrends(keyword = "flu",
                       geo = "US-CA",
                       time = paste(MMWRweek2Date(2006, 41), MMWRweek2Date(2011, 40)))$interest_over_time %>%
  mutate(hits = case_when(hits == "<1" ~ 1,
                          TRUE ~ as.numeric(hits)))

CA_flu_1015 <- gtrends(keyword = "flu",
                       geo = "US-CA",
                       time = paste(MMWRweek2Date(2010, 41), MMWRweek2Date(2015, 40)))$interest_over_time %>%
  mutate(hits = case_when(hits == "<1" ~ 1,
                          TRUE ~ as.numeric(hits)))

CA_flu_1419 <- gtrends(keyword = "flu",
                       geo = "US-CA")$interest_over_time %>%
  mutate(hits = case_when(hits == "<1" ~ 1,
                          TRUE ~ as.numeric(hits)))

# Set up ratios to normalize all Gtrends data to scale from 2010-2015
CA_gratio_0607 <- US_flu_ratio(CA_flu_0407, CA_flu_0611)
CA_gratio_1011 <- US_flu_ratio(CA_flu_0611, CA_flu_1015)
CA_inv_gratio_1415 <- US_flu_ratio(CA_flu_1419, CA_flu_1015)

# Merge Gtrends data and rescale to 2010-2015 scale
gtrend_CA_flu_merge <- filter(CA_flu_0407, !date %in% CA_flu_0611$date) %>%
  mutate(hits = hits / CA_gratio_0607) %>%
  bind_rows(CA_flu_0611) %>%
  filter(!date %in% CA_flu_1015$date) %>%
  mutate(hits = hits / CA_gratio_1011) %>%
  bind_rows(CA_flu_1015) %>%
  filter(!date %in% CA_flu_1419$date) %>%
  bind_rows(CA_flu_1419 %>%
              mutate(hits = hits / CA_inv_gratio_1415)) %>%
  select(date, hits) %>%
  mutate(week = MMWRweek(date)[[2]],
         year = MMWRweek(date)[[1]],
         season = ifelse(week >= 40,
                         paste0(year, "/", year + 1),
                         paste0(year - 1, "/", year)),
         hits = case_when(near(hits, 0) ~ 1,
                          TRUE ~ hits))

# Region 10 - Washington
WA_flu_0407 <- gtrends(keyword = "flu",
                       geo = "US-WA",
                       time = paste(MMWRweek2Date(2004, 41), MMWRweek2Date(2007, 40)))$interest_over_time %>%
  mutate(hits = case_when(hits == "<1" ~ 1,
                          TRUE ~ as.numeric(hits)))

WA_flu_0611 <- gtrends(keyword = "flu",
                       geo = "US-WA",
                       time = paste(MMWRweek2Date(2006, 41), MMWRweek2Date(2011, 40)))$interest_over_time %>%
  mutate(hits = case_when(hits == "<1" ~ 1,
                          TRUE ~ as.numeric(hits)))

WA_flu_1015 <- gtrends(keyword = "flu",
                       geo = "US-WA",
                       time = paste(MMWRweek2Date(2010, 41), MMWRweek2Date(2015, 40)))$interest_over_time %>%
  mutate(hits = case_when(hits == "<1" ~ 1,
                          TRUE ~ as.numeric(hits)))

WA_flu_1419 <- gtrends(keyword = "flu",
                       geo = "US-WA")$interest_over_time %>%
  mutate(hits = case_when(hits == "<1" ~ 1,
                          TRUE ~ as.numeric(hits)))

# Set up ratios to normalize all Gtrends data to scale from 2010-2015
WA_gratio_0607 <- US_flu_ratio(WA_flu_0407, WA_flu_0611)
WA_gratio_1011 <- US_flu_ratio(WA_flu_0611, WA_flu_1015)
WA_inv_gratio_1415 <- US_flu_ratio(WA_flu_1419, WA_flu_1015)

# Merge Gtrends data and rescale to 2010-2015 scale
gtrend_WA_flu_merge <- filter(WA_flu_0407, !date %in% WA_flu_0611$date) %>%
  mutate(hits = hits / WA_gratio_0607) %>%
  bind_rows(WA_flu_0611) %>%
  filter(!date %in% WA_flu_1015$date) %>%
  mutate(hits = hits / WA_gratio_1011) %>%
  bind_rows(WA_flu_1015) %>%
  filter(!date %in% WA_flu_1419$date) %>%
  bind_rows(WA_flu_1419 %>%
              mutate(hits = hits / WA_inv_gratio_1415)) %>%
  select(date, hits) %>%
  mutate(week = MMWRweek(date)[[2]],
         year = MMWRweek(date)[[1]],
         season = ifelse(week >= 40,
                         paste0(year, "/", year + 1),
                         paste0(year - 1, "/", year)),
         hits = case_when(near(hits, 0) ~ 1,
                          TRUE ~ hits))

save(gtrend_MA_flu_merge, gtrend_NY_flu_merge, gtrend_PA_flu_merge,
     gtrend_FL_flu_merge, gtrend_IL_flu_merge, gtrend_TX_flu_merge,
     gtrend_MO_flu_merge, gtrend_CO_flu_merge, gtrend_CA_flu_merge,
     gtrend_WA_flu_merge, file = "Data/state_gtrend.Rdata")

## Explore other search topics based on Google Correlate -------
# 'get over the flu'
getoverflu_1014<- gtrends(keyword = "get over the flu",
                      geo = "US",
                      time = paste(MMWRweek2Date(2010, 40), MMWRweek2Date(2014, 40)))$interest_over_time

getoverflu_1418 <- gtrends(keyword = "get over the flu",
                       geo = c("US"))$interest_over_time

getoverflu_merge <- inner_join(getoverflu_1014 %>%
                             select(date, old_hits = hits),
                           getoverflu_1418 %>%
                             select(date, new_hits = hits),
                           by = "date") %>%
  mutate(ratio = old_hits / new_hits)

ggplot(getoverflu_merge) +
  geom_line(aes(x = date, y = old_hits), color = "red") +
  geom_line(aes(x = date, y = new_hits), color = "blue")

ggplot(getoverflu_merge) +
  geom_line(aes(x = date, y = ratio))

# 'influenza type a'
flutypea_1014<- gtrends(keyword = "influenza type a",
                          geo = "US",
                          time = paste(MMWRweek2Date(2010, 40), MMWRweek2Date(2014, 40)))$interest_over_time

flutypea_1418 <- gtrends(keyword = "influenza type a",
                           geo = c("US"))$interest_over_time

flutypea_merge <- inner_join(flutypea_1014 %>%
                               select(date, old_hits = hits) %>%
                               filter(date %within% interval(MMWRweek2Date(2013, 42), 
                                                             MMWRweek2Date(2014, 18))),
                             flutypea_1418 %>%
                               select(date, new_hits = hits) %>%
                               filter(date %within% interval(MMWRweek2Date(2013, 42), 
                                                             MMWRweek2Date(2014, 18))),
                               by = "date") %>%
  mutate(ratio = old_hits / new_hits)
?lubridate::interval
ggplot(flutypea_merge) +
  geom_line(aes(x = date, y = old_hits), color = "red") +
  geom_line(aes(x = date, y = new_hits), color = "blue")

ggplot(flutypea_merge) +
  geom_line(aes(x = date, y = ratio))

