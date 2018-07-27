require(tidyverse)
require(gtrendsR)
require(lubridate)
require(MMWRweek)

# Create conversion ratios for Google Trends in each region
gtrend_ratio <- function(location) {
  flu_1014<- gtrends(keyword = "influenza",
                     geo = location,
                     time = paste(MMWRweek2Date(2010, 40), MMWRweek2Date(2014, 40)))$interest_over_time

  flu_1418 <- gtrends(keyword = "influenza",
                      geo = location)$interest_over_time
  
  flu_merge <- inner_join(flu_1014 %>%
                             select(date, old_hits = hits) %>%
                             filter(date %within% interval(MMWRweek2Date(2013, 42), 
                                                           MMWRweek2Date(2014, 18))),
                           flu_1418 %>%
                             select(date, new_hits = hits) %>%
                             filter(date %within% interval(MMWRweek2Date(2013, 42), 
                                                           MMWRweek2Date(2014, 18))),
                           by = "date") %>%
  mutate(logratio = log2(old_hits / new_hits)) 
  
  2^mean(flu_merge$logratio)
}
