# Google Trends data
library(tidyverse)
library(MMWRweek)
library(gtrendsR)
library(lubridate)

source("R/utils.R")

# Fetch Google Trends data -----
gtrend_US_flu_merge <- fetch_gtrend("US")

save(gtrend_US_flu_merge, file = "Data/Gtrends.Rdata")

# Region 1 - Massachusetts
gtrend_MA_flu_merge <- fetch_gtrend("MA")
# Region 2 - New York
gtrend_NY_flu_merge <- fetch_gtrend("NY")
# Region 3 - Pennsylvania
gtrend_PA_flu_merge <- fetch_gtrend("PA")
gtrend_DE_flu_merge <- fetch_gtrend("DE")
gtrend_MD_flu_merge <- fetch_gtrend("MD")
gtrend_VA_flu_merge <- fetch_gtrend("VA")
gtrend_WV_flu_merge <- fetch_gtrend("WV")
# Region 4 - Florida
gtrend_FL_flu_merge <- fetch_gtrend("FL")
# Region 5 - Illinois
gtrend_IL_flu_merge <- fetch_gtrend("IL")
# Region 6 - Texas
gtrend_TX_flu_merge <- fetch_gtrend("TX")
# Region 7 - Missouri
gtrend_MO_flu_merge <- fetch_gtrend("MO")
# Region 8 - Colorado
gtrend_CO_flu_merge <- fetch_gtrend("CO")
# Region 9 - California
gtrend_CA_flu_merge <- fetch_gtrend("CA")
# Region 10 - Washington
gtrend_WA_flu_merge <- fetch_gtrend("WA")

save(gtrend_MA_flu_merge, gtrend_NY_flu_merge, gtrend_PA_flu_merge,
     gtrend_FL_flu_merge, gtrend_IL_flu_merge, gtrend_TX_flu_merge,
     gtrend_MO_flu_merge, gtrend_CO_flu_merge, gtrend_CA_flu_merge,
     gtrend_WA_flu_merge, gtrend_DE_flu_merge, gtrend_MD_flu_merge,
     gtrend_VA_flu_merge, gtrend_WV_flu_merge,
     file = "Data/state_gtrend.Rdata")

