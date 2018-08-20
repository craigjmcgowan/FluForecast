# Generate scores
library(tidyverse)
library(FluSight)

source("R/read_forecasts.R")
# Read in retrospective forecasts -------
subtype_forecasts <- read_forecasts("Forecasts/Subtype Historical Average")
verify_entry_file("Forecasts/Subtype Historical Average/2012-2013/EW1.csv")

dir <- "Forecasts/Subtype Historical Average"
