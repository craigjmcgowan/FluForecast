# Google Trends data
library(tidyverse)
library(MMWRweek)
library(gtrendsR)
library(lubridate)

source("R/utils.R")

state_matchup <- tibble(state_abb = state.abb,
                        location = state.name) 

# Fetch Google Trends data -----
gtrend_US_flu_merge <- fetch_gtrend("US")

save(gtrend_US_flu_merge, file = "Data/Gtrends.Rdata")

# State Google Trends data
gtrend_state_list <- tibble(state_abb = state.abb) %>%
  mutate(data = map(state_abb, ~ fetch_gtrend(.)))

gtrend_state_flu_merge <- unnest(gtrend_state_list) %>%
  full_join(state_matchup, by = "state_abb") %>%
  select(-state_abb)

# # Region 1 - Massachusetts
# gtrend_MA_flu_merge <- filter(gtrend_state_list, state_abb == "MA") %>%
#   select(-state_abb)
# gtrend_CT_flu_merge <- filter(gtrend_state_list, state_abb == "CT") %>%
#   select(-state_abb)
# gtrend_RI_flu_merge <- filter(gtrend_state_list, state_abb == "RI") %>%
#   select(-state_abb)
# gtrend_VT_flu_merge <- filter(gtrend_state_list, state_abb == "VT") %>%
#   select(-state_abb)
# gtrend_NH_flu_merge <- filter(gtrend_state_list, state_abb == "NH") %>%
#   select(-state_abb)
# gtrend_ME_flu_merge <- filter(gtrend_state_list, state_abb == "ME") %>%
#   select(-state_abb)
# 
# # Region 2 - New York
# gtrend_NY_flu_merge <- filter(gtrend_state_list, state_abb == "NY") %>%
#   select(-state_abb)
# gtrend_NJ_flu_merge <- filter(gtrend_state_list, state_abb == "NJ") %>%
#   select(-state_abb)
# 
# # Region 3 - Pennsylvania
# gtrend_PA_flu_merge <- filter(gtrend_state_list, state_abb == "PA") %>%
#   select(-state_abb)
# gtrend_DE_flu_merge <- filter(gtrend_state_list, state_abb == "DE") %>%
#   select(-state_abb)
# gtrend_MD_flu_merge <- filter(gtrend_state_list, state_abb == "MD") %>%
#   select(-state_abb)
# gtrend_VA_flu_merge <- filter(gtrend_state_list, state_abb == "VA") %>%
#   select(-state_abb)
# gtrend_WV_flu_merge <- filter(gtrend_state_list, state_abb == "WV") %>%
#   select(-state_abb)
# 
# # Region 4 - Florida
# gtrend_NC_flu_merge <- filter(gtrend_state_list, state_abb == "NC") %>%
#   select(-state_abb)
# gtrend_SC_flu_merge <- filter(gtrend_state_list, state_abb == "SC") %>%
#   select(-state_abb)
# gtrend_GA_flu_merge <- filter(gtrend_state_list, state_abb == "GA") %>%
#   select(-state_abb)
# gtrend_FL_flu_merge <- filter(gtrend_state_list, state_abb == "FL") %>%
#   select(-state_abb)
# gtrend_AL_flu_merge <- filter(gtrend_state_list, state_abb == "AL") %>%
#   select(-state_abb)
# gtrend_MS_flu_merge <- filter(gtrend_state_list, state_abb == "MS") %>%
#   select(-state_abb)
# gtrend_KY_flu_merge <- filter(gtrend_state_list, state_abb == "KY") %>%
#   select(-state_abb)
# gtrend_TN_flu_merge <- filter(gtrend_state_list, state_abb == "TN") %>%
#   select(-state_abb)
# 
# # Region 5 - Illinois
# gtrend_IN_flu_merge <- filter(gtrend_state_list, state_abb == "IN") %>%
#   select(-state_abb)
# gtrend_IL_flu_merge <- filter(gtrend_state_list, state_abb == "IL") %>%
#   select(-state_abb)
# gtrend_OH_flu_merge <- filter(gtrend_state_list, state_abb == "OH") %>%
#   select(-state_abb)
# gtrend_MI_flu_merge <- filter(gtrend_state_list, state_abb == "MI") %>%
#   select(-state_abb)
# gtrend_WI_flu_merge <- filter(gtrend_state_list, state_abb == "WI") %>%
#   select(-state_abb)
# gtrend_MN_flu_merge <- filter(gtrend_state_list, state_abb == "MN") %>%
#   select(-state_abb)
# 
# # Region 6 - Texas
# gtrend_LA_flu_merge <- filter(gtrend_state_list, state_abb == "LA") %>%
#   select(-state_abb)
# gtrend_AR_flu_merge <- filter(gtrend_state_list, state_abb == "AR") %>%
#   select(-state_abb)
# gtrend_OK_flu_merge <- filter(gtrend_state_list, state_abb == "OK") %>%
#   select(-state_abb)
# gtrend_TX_flu_merge <- filter(gtrend_state_list, state_abb == "TX") %>%
#   select(-state_abb)
# gtrend_NM_flu_merge <- filter(gtrend_state_list, state_abb == "NM") %>%
#   select(-state_abb)
# 
# # Region 7 - Missouri
# gtrend_IA_flu_merge <- filter(gtrend_state_list, state_abb == "IA") %>%
#   select(-state_abb)
# gtrend_MO_flu_merge <- filter(gtrend_state_list, state_abb == "MO") %>%
#   select(-state_abb)
# gtrend_NE_flu_merge <- filter(gtrend_state_list, state_abb == "NE") %>%
#   select(-state_abb)
# gtrend_KS_flu_merge <- filter(gtrend_state_list, state_abb == "KS") %>%
#   select(-state_abb)
# 
# # Region 8 - Colorado
# gtrend_ND_flu_merge <- filter(gtrend_state_list, state_abb == "ND") %>%
#   select(-state_abb)
# gtrend_SD_flu_merge <- filter(gtrend_state_list, state_abb == "SD") %>%
#   select(-state_abb)
# gtrend_MT_flu_merge <- filter(gtrend_state_list, state_abb == "MT") %>%
#   select(-state_abb) 
# gtrend_WY_flu_merge <- filter(gtrend_state_list, state_abb == "WY") %>%
#   select(-state_abb)
# gtrend_CO_flu_merge <- filter(gtrend_state_list, state_abb == "CO") %>%
#   select(-state_abb)
# gtrend_UT_flu_merge <- filter(gtrend_state_list, state_abb == "UT") %>%
#   select(-state_abb)
# 
# # Region 9 - California
# gtrend_CA_flu_merge <- filter(gtrend_state_list, state_abb == "CA") %>%
#   select(-state_abb)
# gtrend_NV_flu_merge <- filter(gtrend_state_list, state_abb == "NV") %>%
#   select(-state_abb)
# gtrend_AZ_flu_merge <- filter(gtrend_state_list, state_abb == "AZ") %>%
#   select(-state_abb)
# gtrend_HI_flu_merge <- filter(gtrend_state_list, state_abb == "HI") %>%
#   select(-state_abb)
# 
# # Region 10 - Washington
# gtrend_WA_flu_merge <- filter(gtrend_state_list, state_abb == "WA") %>%
#   select(-state_abb)
# gtrend_ID_flu_merge <- filter(gtrend_state_list, state_abb == "ID") %>%
#   select(-state_abb)
# gtrend_OR_flu_merge <- filter(gtrend_state_list, state_abb == "OR") %>%
#   select(-state_abb)
# gtrend_AK_flu_merge <- filter(gtrend_state_list, state_abb == "AK") %>%
#   select(-state_abb)

save(gtrend_state_flu_merge,
     file = "Data/state_gtrend.Rdata")

