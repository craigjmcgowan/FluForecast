# Google Trends data
library(tidyverse)
library(MMWRweek)
library(gtrendsR)


# Try overall US
US_flu_1014<- gtrends(keyword = "influenza",
                      geo = "US",
                      time = paste(MMWRweek2Date(2010, 40), MMWRweek2Date(2014, 40)))$interest_over_time

US_flu_1418 <- gtrends(keyword = "influenza",
                       geo = "US")$interest_over_time

US_flu_merge <- inner_join(US_flu_1014 %>%
                             select(date, old_hits = hits),
                           US_flu_1418 %>%
                             select(date, new_hits = hits),
                           by = "date") %>%
  mutate(ratio = old_hits / new_hits) 

ggplot(US_flu_merge) +
  geom_line(aes(x = date, y = old_hits), color = "red") +
  geom_line(aes(x = date, y = new_hits), color = "blue")

ggplot(US_flu_merge, aes(x = date, y = ratio)) +
  geom_line()


# Explore trends by state - 1 state in each region ------------

# Region 1 - Massachusetts
MA_flu_1014<- gtrends(keyword = "influenza",
                      geo = "US-MA",
                      time = paste(MMWRweek2Date(2010, 40), MMWRweek2Date(2014, 40)))$interest_over_time

MA_flu_1418 <- gtrends(keyword = "influenza",
                       geo = "US-MA")$interest_over_time

MA_flu_merge <- full_join(MA_flu_1014 %>%
                            select(date, old_hits = hits),
                          MA_flu_1418 %>%
                            select(date, new_hits = hits),
                          by = "date") %>%
  mutate(ratio = old_hits / new_hits) %>%
  filter(!is.na(ratio))

ggplot(MA_flu_merge) +
  geom_line(aes(x = date, y = old_hits), color = "red") +
  geom_line(aes(x = date, y = new_hits), color = "blue")
# Not too bad - matches up pretty well

ggplot(MA_flu_merge, aes(x = date, y = ratio)) +
  geom_line()

# Region 2 - New York
NY_flu_1014<- gtrends(keyword = "influenza",
                      geo = "US-NY",
                      time = paste(MMWRweek2Date(2010, 40), MMWRweek2Date(2014, 40)))$interest_over_time

NY_flu_1418 <- gtrends(keyword = "influenza",
                       geo = "US-NY")$interest_over_time

NY_flu_merge <- inner_join(NY_flu_1014 %>%
                            select(date, old_hits = hits),
                          NY_flu_1418 %>%
                            select(date, new_hits = hits),
                          by = "date") %>%
  mutate(ratio = old_hits / new_hits)

ggplot(NY_flu_merge) +
  geom_line(aes(x = date, y = old_hits), color = "red") +
  geom_line(aes(x = date, y = new_hits), color = "blue")
# OK, not great

ggplot(NY_flu_merge, aes(x = date, y = ratio)) +
  geom_line()

# Region 3 - Pennsylvania
PA_flu_1014<- gtrends(keyword = "influenza",
                      geo = "US-PA",
                      time = paste(MMWRweek2Date(2010, 40), MMWRweek2Date(2014, 40)))$interest_over_time

PA_flu_1418 <- gtrends(keyword = "influenza",
                       geo = "US-PA")$interest_over_time

PA_flu_merge <- inner_join(PA_flu_1014 %>%
                             select(date, old_hits = hits),
                           PA_flu_1418 %>%
                             select(date, new_hits = hits),
                           by = "date") %>%
  mutate(ratio = old_hits / new_hits)

ggplot(PA_flu_merge) +
  geom_line(aes(x = date, y = old_hits), color = "red") +
  geom_line(aes(x = date, y = new_hits), color = "blue")
# OK, not great

ggplot(PA_flu_merge, aes(x = date, y = ratio)) +
  geom_line()

# Region 4 - Florida
FL_flu_1014<- gtrends(keyword = "influenza",
                      geo = "US-FL",
                      time = paste(MMWRweek2Date(2010, 40), MMWRweek2Date(2014, 40)))$interest_over_time

FL_flu_1418 <- gtrends(keyword = "influenza",
                       geo = "US-FL")$interest_over_time

FL_flu_merge <- inner_join(FL_flu_1014 %>%
                             select(date, old_hits = hits),
                           FL_flu_1418 %>%
                             select(date, new_hits = hits),
                           by = "date") %>%
  mutate(ratio = old_hits / new_hits)

ggplot(FL_flu_merge) +
  geom_line(aes(x = date, y = old_hits), color = "red") +
  geom_line(aes(x = date, y = new_hits), color = "blue")
# Pretty good pattern match

ggplot(FL_flu_merge, aes(x = date, y = ratio)) +
  geom_line()

# Region 5 - Illinois
IL_flu_1014<- gtrends(keyword = "influenza",
                      geo = "US-IL",
                      time = paste(MMWRweek2Date(2010, 40), MMWRweek2Date(2014, 40)))$interest_over_time

IL_flu_1418 <- gtrends(keyword = "influenza",
                       geo = "US-IL")$interest_over_time

IL_flu_merge <- inner_join(IL_flu_1014 %>%
                             select(date, old_hits = hits),
                           IL_flu_1418 %>%
                             select(date, new_hits = hits),
                           by = "date") %>%
  mutate(ratio = old_hits / new_hits)

ggplot(IL_flu_merge) +
  geom_line(aes(x = date, y = old_hits), color = "red") +
  geom_line(aes(x = date, y = new_hits), color = "blue")
# OK, not great

ggplot(IL_flu_merge, aes(x = date, y = ratio)) +
  geom_line()

# Region 6 - Texas
TX_flu_1014<- gtrends(keyword = "influenza",
                      geo = "US-TX",
                      time = paste(MMWRweek2Date(2010, 40), MMWRweek2Date(2014, 40)))$interest_over_time

TX_flu_1418 <- gtrends(keyword = "influenza",
                       geo = "US-TX")$interest_over_time

TX_flu_merge <- inner_join(TX_flu_1014 %>%
                             select(date, old_hits = hits),
                           TX_flu_1418 %>%
                             select(date, new_hits = hits),
                           by = "date") %>%
  mutate(ratio = old_hits / new_hits)

ggplot(TX_flu_merge) +
  geom_line(aes(x = date, y = old_hits), color = "red") +
  geom_line(aes(x = date, y = new_hits), color = "blue")
# Pretty good match up

ggplot(TX_flu_merge, aes(x = date, y = ratio)) +
  geom_line()

# Region 7 - Missouri
MO_flu_1014<- gtrends(keyword = "influenza",
                      geo = "US-MO",
                      time = paste(MMWRweek2Date(2010, 40), MMWRweek2Date(2014, 40)))$interest_over_time

MO_flu_1418 <- gtrends(keyword = "influenza",
                       geo = "US-MO")$interest_over_time

MO_flu_merge <- inner_join(MO_flu_1014 %>%
                             select(date, old_hits = hits),
                           MO_flu_1418 %>%
                             select(date, new_hits = hits),
                           by = "date") %>%
  mutate(ratio = old_hits / new_hits)

ggplot(MO_flu_merge) +
  geom_line(aes(x = date, y = old_hits), color = "red") +
  geom_line(aes(x = date, y = new_hits), color = "blue")
# Pretty good match up

ggplot(MO_flu_merge, aes(x = date, y = ratio)) +
  geom_line()

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
