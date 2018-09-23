# Google Trends data
library(tidyverse)
library(MMWRweek)
library(gtrendsR)
library(lubridate)

<<<<<<< HEAD
load("Data/Gtrends.Rdata")
load("Data/virologic.Rdata")
load("Data/ili.Rdata")
source("R/utils.R")

# Merge Google Trends and ILI
gtrend_ILI_merge <- full_join(
  gtrend_US_flu_merge, 
  filter(ili_current, location == "US National") %>%
    select(season, week, year, ILI),
  by = c("season", "week", "year")
  ) %>% 
  full_join(filter(virologic_combined, location == "US National") %>%
    select(season, week, year, h1_per_samples, h3_per_samples, b_per_samples),
    by = c("season", "week", "year")) %>%
  filter(year >= 2004, season != "2003/2004")

# Plot of US Google Trend data since 2004
ggplot(mutate(gtrend_US_flu_merge, order_week = week_inorder(week, season)) %>%
         filter(order_week <= 72)) +
  geom_line(aes(x = order_week, y = hits, color = season)) +
  scale_y_continuous(limits = c(0, 100))

# Convert to ts() object
gtrend_US_ts <- ts(gtrend_US_flu_merge$hits, frequency = 52, start = c(2004, 40))
=======

# Try overall US
US_flu_1014<- gtrends(keyword = "influenza",
                      geo = "US",
                      time = paste(MMWRweek2Date(2010, 40), MMWRweek2Date(2014, 40)))$interest_over_time

US_flu_1418 <- gtrends(keyword = "influenza",
                       geo = "US")$interest_over_time

US_flu_merge <- inner_join(US_flu_1014 %>%
                             select(date, old_hits = hits) %>%
                             filter(date %within% interval(MMWRweek2Date(2013, 42), 
                                                           MMWRweek2Date(2014, 18))),
                           US_flu_1418 %>%
                             select(date, new_hits = hits) %>%
                             filter(date %within% interval(MMWRweek2Date(2013, 42), 
                                                           MMWRweek2Date(2014, 18))),
                           by = "date") %>%
  mutate(ratio = old_hits / new_hits,
         logratio = log2(ratio)) 

ggplot(US_flu_merge) +
  geom_line(aes(x = date, y = old_hits), color = "red") +
  geom_line(aes(x = date, y = new_hits), color = "blue")
>>>>>>> c4102c6aac3ee17dae5abfb2ad68ba2cee545e53

ggplot(US_flu_merge) +
  geom_line(aes(x = date, y = ratio), color = "black") +
  geom_line(aes(x = date, y = logratio), color = "red")


# Explore trends by state - 1 state in each region ------------

# Region 1 - Massachusetts
<<<<<<< HEAD
MA_flu_1014<- gtrends(keyword = "flu",
                      geo = "US-MA",
                      time = paste(MMWRweek2Date(2010, 40), MMWRweek2Date(2014, 40)))$interest_over_time

MA_flu_1418 <- gtrends(keyword = "flu",
=======
MA_flu_1014<- gtrends(keyword = "influenza",
                      geo = "US-MA",
                      time = paste(MMWRweek2Date(2010, 40), MMWRweek2Date(2014, 40)))$interest_over_time

MA_flu_1418 <- gtrends(keyword = "influenza",
>>>>>>> c4102c6aac3ee17dae5abfb2ad68ba2cee545e53
                       geo = "US-MA")$interest_over_time

MA_flu_merge <- full_join(MA_flu_1014 %>%
                            select(date, old_hits = hits) %>%
                            filter(date %within% interval(MMWRweek2Date(2013, 42), 
                                                          MMWRweek2Date(2014, 18))),
                          MA_flu_1418 %>%
                            select(date, new_hits = hits) %>%
                            filter(date %within% interval(MMWRweek2Date(2013, 42), 
                                                          MMWRweek2Date(2014, 18))),
                          by = "date") %>%
  mutate(ratio = old_hits / new_hits,
         logratio = log2(ratio)) %>%
  filter(!is.na(ratio))

<<<<<<< HEAD
ggplot(MA_flu_1418) +
  geom_line(aes(x = date, y = hits))

=======
>>>>>>> c4102c6aac3ee17dae5abfb2ad68ba2cee545e53
ggplot(MA_flu_merge) +
  geom_line(aes(x = date, y = old_hits), color = "red") +
  geom_line(aes(x = date, y = new_hits), color = "blue")
# Not too bad - matches up pretty well

ggplot(MA_flu_merge) +
  geom_line(aes(x = date, y = ratio), color = "black") +
  geom_line(aes(x = date, y = logratio), color = "red")

# Region 2 - New York
NY_flu_1014<- gtrends(keyword = "influenza",
                      geo = "US-NY",
                      time = paste(MMWRweek2Date(2010, 40), MMWRweek2Date(2014, 40)))$interest_over_time

NY_flu_1418 <- gtrends(keyword = "influenza",
                       geo = "US-NY")$interest_over_time

NY_flu_merge <- inner_join(NY_flu_1014 %>%
                            select(date, old_hits = hits) %>%
                             filter(date %within% interval(MMWRweek2Date(2013, 42), 
                                                           MMWRweek2Date(2014, 18))),
                          NY_flu_1418 %>%
                            select(date, new_hits = hits) %>%
                            filter(date %within% interval(MMWRweek2Date(2013, 42), 
                                                          MMWRweek2Date(2014, 18))),
                          by = "date") %>%
  mutate(ratio = old_hits / new_hits,
         logratio = log2(ratio))

ggplot(NY_flu_merge) +
  geom_line(aes(x = date, y = old_hits), color = "red") +
  geom_line(aes(x = date, y = new_hits), color = "blue")
# OK, not great

ggplot(NY_flu_merge) +
  geom_line(aes(x = date, y = ratio), color = "black") +
  geom_line(aes(x = date, y = logratio), color = "red")

# Region 3 - Pennsylvania
PA_flu_1014<- gtrends(keyword = "influenza",
                      geo = "US-PA",
                      time = paste(MMWRweek2Date(2010, 40), MMWRweek2Date(2014, 40)))$interest_over_time

PA_flu_1418 <- gtrends(keyword = "influenza",
                       geo = "US-PA")$interest_over_time

PA_flu_merge <- inner_join(PA_flu_1014 %>%
                             select(date, old_hits = hits) %>%
                             filter(date %within% interval(MMWRweek2Date(2013, 42), 
                                                           MMWRweek2Date(2014, 18))),
                           PA_flu_1418 %>%
                             select(date, new_hits = hits) %>%
                             filter(date %within% interval(MMWRweek2Date(2013, 42), 
                                                           MMWRweek2Date(2014, 18))),
                           by = "date") %>%
  mutate(ratio = old_hits / new_hits,
         logratio = log2(ratio))

ggplot(PA_flu_merge) +
  geom_line(aes(x = date, y = old_hits), color = "red") +
  geom_line(aes(x = date, y = new_hits), color = "blue")
# OK, not great

ggplot(PA_flu_merge) +
  geom_line(aes(x = date, y = ratio), color = "black") +
  geom_line(aes(x = date, y = logratio), color = "red")

# Region 4 - Florida
FL_flu_1014<- gtrends(keyword = "influenza",
                      geo = "US-FL",
                      time = paste(MMWRweek2Date(2010, 40), MMWRweek2Date(2014, 40)))$interest_over_time

FL_flu_1418 <- gtrends(keyword = "influenza",
                       geo = "US-FL")$interest_over_time

FL_flu_merge <- inner_join(FL_flu_1014 %>%
                             select(date, old_hits = hits) %>%
                             filter(date %within% interval(MMWRweek2Date(2013, 42), 
                                                           MMWRweek2Date(2014, 18))),
                           FL_flu_1418 %>%
                             select(date, new_hits = hits) %>%
                             filter(date %within% interval(MMWRweek2Date(2013, 42), 
                                                           MMWRweek2Date(2014, 18))),
                           by = "date") %>%
  mutate(ratio = old_hits / new_hits,
         logratio = log2(ratio))

ggplot(FL_flu_merge) +
  geom_line(aes(x = date, y = old_hits), color = "red") +
  geom_line(aes(x = date, y = new_hits), color = "blue")
# Pretty good pattern match

ggplot(FL_flu_merge) +
  geom_line(aes(x = date, y = ratio), color = "black") +
  geom_line(aes(x = date, y = logratio), color = "red")

# Region 5 - Illinois
IL_flu_1014<- gtrends(keyword = "influenza",
                      geo = "US-IL",
                      time = paste(MMWRweek2Date(2010, 40), MMWRweek2Date(2014, 40)))$interest_over_time

IL_flu_1418 <- gtrends(keyword = "influenza",
                       geo = "US-IL")$interest_over_time

IL_flu_merge <- inner_join(IL_flu_1014 %>%
                             select(date, old_hits = hits) %>%
                             filter(date %within% interval(MMWRweek2Date(2013, 42), 
                                                           MMWRweek2Date(2014, 18))),
                           IL_flu_1418 %>%
                             select(date, new_hits = hits) %>%
                             mutate(new_hits = ifelse(new_hits == "<1", as.numeric(0),
                                                  as.numeric(new_hits))) %>%
                             filter(date %within% interval(MMWRweek2Date(2013, 42), 
                                                           MMWRweek2Date(2014, 18))),
                           by = "date") %>%
  mutate(ratio = old_hits / new_hits,
         logratio = log2(ratio))

ggplot(IL_flu_merge) +
  geom_line(aes(x = date, y = old_hits), color = "red") +
  geom_line(aes(x = date, y = new_hits), color = "blue")
# OK, not great

ggplot(IL_flu_merge) +
  geom_line(aes(x = date, y = ratio), color = "black") +
  geom_line(aes(x = date, y = logratio), color = "red")

# Region 6 - Texas
TX_flu_1014<- gtrends(keyword = "influenza",
                      geo = "US-TX",
                      time = paste(MMWRweek2Date(2010, 40), MMWRweek2Date(2014, 40)))$interest_over_time

TX_flu_1418 <- gtrends(keyword = "influenza",
                       geo = "US-TX")$interest_over_time

TX_flu_merge <- inner_join(TX_flu_1014 %>%
                             select(date, old_hits = hits) %>%
                             filter(date %within% interval(MMWRweek2Date(2013, 42), 
                                                           MMWRweek2Date(2014, 18))),
                           TX_flu_1418 %>%
                             select(date, new_hits = hits) %>%
                             filter(date %within% interval(MMWRweek2Date(2013, 42), 
                                                           MMWRweek2Date(2014, 18))),
                           by = "date") %>%
  mutate(ratio = old_hits / new_hits,
         logratio = log2(ratio))

ggplot(TX_flu_merge) +
  geom_line(aes(x = date, y = old_hits), color = "red") +
  geom_line(aes(x = date, y = new_hits), color = "blue")
# Pretty good match up

ggplot(TX_flu_merge) +
  geom_line(aes(x = date, y = ratio), color = "black") +
  geom_line(aes(x = date, y = logratio), color = "red")

# Region 7 - Missouri
MO_flu_1014<- gtrends(keyword = "influenza",
                      geo = "US-MO",
                      time = paste(MMWRweek2Date(2010, 40), MMWRweek2Date(2014, 40)))$interest_over_time

MO_flu_1418 <- gtrends(keyword = "influenza",
                       geo = "US-MO")$interest_over_time

MO_flu_merge <- inner_join(MO_flu_1014 %>%
                             select(date, old_hits = hits) %>%
                             filter(date %within% interval(MMWRweek2Date(2013, 42), 
                                                           MMWRweek2Date(2014, 18))),
                           MO_flu_1418 %>%
                             select(date, new_hits = hits) %>%
                             filter(date %within% interval(MMWRweek2Date(2013, 42), 
                                                           MMWRweek2Date(2014, 18))),
                           by = "date") %>%
  mutate(ratio = old_hits / new_hits,
         logratio = log2(ratio))

ggplot(MO_flu_merge) +
  geom_line(aes(x = date, y = old_hits), color = "red") +
  geom_line(aes(x = date, y = new_hits), color = "blue")
# Pretty good match up

ggplot(MO_flu_merge) +
  geom_line(aes(x = date, y = ratio), color = "black") +
  geom_line(aes(x = date, y = logratio), color = "red")

# Region 8 - Colorado
CO_flu_1014<- gtrends(keyword = "influenza",
                      geo = "US-CO",
                      time = paste(MMWRweek2Date(2010, 40), MMWRweek2Date(2014, 40)))$interest_over_time

CO_flu_1418 <- gtrends(keyword = "influenza",
                       geo = "US-CO")$interest_over_time

CO_flu_merge <- inner_join(CO_flu_1014 %>%
                             select(date, old_hits = hits),
                           CO_flu_1418 %>%
                             select(date, new_hits = hits),
                           by = "date") %>%
  mutate(ratio = old_hits / new_hits)

ggplot(CO_flu_merge) +
  geom_line(aes(x = date, y = old_hits), color = "red") +
  geom_line(aes(x = date, y = new_hits), color = "blue")
# Pretty good match up

ggplot(CO_flu_merge, aes(x = date, y = ratio)) +
  geom_line()

# Region 9 - California
CA_flu_1014<- gtrends(keyword = "influenza",
                      geo = "US-CA",
                      time = paste(MMWRweek2Date(2010, 40), MMWRweek2Date(2014, 40)))$interest_over_time

CA_flu_1418 <- gtrends(keyword = "influenza",
                       geo = "US-CA")$interest_over_time

CA_flu_merge <- inner_join(CA_flu_1014 %>%
                             select(date, old_hits = hits),
                           CA_flu_1418 %>%
                             select(date, new_hits = hits),
                           by = "date") %>%
  mutate(ratio = old_hits / new_hits)

ggplot(CA_flu_merge) +
  geom_line(aes(x = date, y = old_hits), color = "red") +
  geom_line(aes(x = date, y = new_hits), color = "blue")
# Pretty good match up

ggplot(CA_flu_merge, aes(x = date, y = ratio)) +
  geom_line()

# Region 10 - Washington
WA_flu_1014<- gtrends(keyword = "influenza",
                      geo = "US-WA",
                      time = paste(MMWRweek2Date(2010, 40), MMWRweek2Date(2014, 40)))$interest_over_time

WA_flu_1418 <- gtrends(keyword = "influenza",
                       geo = c("US-WA", "US-OR"))$interest_over_time

WA_flu_merge <- inner_join(WA_flu_1014 %>%
                             select(date, old_hits = hits),
                           WA_flu_1418 %>%
                             select(date, new_hits = hits),
                           by = "date") %>%
  mutate(ratio = old_hits / new_hits)

ggplot(WA_flu_merge) +
  geom_line(aes(x = date, y = old_hits), color = "red") +
  geom_line(aes(x = date, y = new_hits), color = "blue")
# Pretty good match up

ggplot(WA_flu_merge, aes(x = date, y = ratio)) +
  geom_line()


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

