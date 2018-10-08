# Log virus subtype data throughout season
library(cdcfluview)

WHO_national <- who_nrevss(region = "national", year = 2018)
WHO_regional <- who_nrevss(region = "hhs", year = 2018)

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
          
          