# Log virus subtype data throughout season
library(cdcfluview)
library(MMWRweek)

# Assign year based on current week
MMWR_week <- MMWRweek(Sys.Date())

year <- dplyr::case_when(
  MMWR_week$MMWRweek[1] < 40 ~ MMWR_week$MMWRyear[1] - 1,
  TRUE ~ MMWR_week$MMWRyear[1]
)

WHO_national <- who_nrevss(region = "national", year = year)
WHO_regional <- who_nrevss(region = "hhs", year = year)
WHO_state <- who_nrevss(region = "state", year = year)

# Save most recent week published
last_week <- str_pad(tail(WHO_national$public_health_labs$week, 1), 2, "left", "0")

# Save data
write_csv(WHO_national$public_health_labs,
          path = paste0("Data/WHO-NREVSS Data/EW", last_week,
                        "Nat_PH_labs.csv"))

write_csv(WHO_national$clinical_labs,
          path = paste0("Data/WHO-NREVSS Data/EW", last_week,
                        "Nat_clin_labs.csv"))

write_csv(WHO_regional$public_health_labs,
          path = paste0("Data/WHO-NREVSS Data/EW", last_week,
                        "Reg_PH_labs.csv"))

write_csv(WHO_regional$clinical_labs,
          path = paste0("Data/WHO-NREVSS Data/EW", last_week,
                        "Reg_clin_labs.csv"))
          
write_csv(WHO_state$public_health_labs,
          path = paste0("Data/WHO-NREVSS Data/EW", last_week,
                        "State_PH_labs.csv"))

write_csv(WHO_state$clinical_labs,
          path = paste0("Data/WHO-NREVSS Data/EW", last_week,
                        "State_clin_labs.csv"))
          