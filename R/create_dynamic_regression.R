# Dynamic regression model 

library(tidyverse)
library(forecast)
library(lubridate)
library(FluSight)
library(MMWRweek)

# Load functions
source("R/utils.R")

# Load data
load("Data/ili.Rdata")
load("Data/virologic.Rdata")
load("Data/Gtrends.Rdata")

# Combine datasets together
flu_data_merge <- select(ili_current, epiweek, ILI, year, week, season, location) %>%
  inner_join(select(virologic_combined, location, season, year, week, h1_per_samples,
                   h3_per_samples, b_per_samples),
            by = c("location", "season", "year", "week")) %>%
  inner_join(gtrend_US_flu_merge,
            by = c("season", "year", "week")) %>%
  # Remove 2008/2009 and 2009/2010 seasons due to pandemic activity
  filter(!season %in% c("2008/2009", "2009/2010")) %>%
  # Remove week 33 in 2014 so all seasons have 52 weeks - minimal activity
  filter(!(year == 2014 & week == 33)) %>%
  # Keep only US for now
  filter(location == "US National")


### Decisions to make ###

## Number of Fourier terms
## Structure of ARIMA errors
## What virologic/Gtrend data to include
## If/how to include estimates of backfill

## Evaluate based on CDC log score


## Plan of action:

## Set up model of some kind to produce CDC-style forecast sequentially for each week of the season
## Set up way to score it
## Scale up to evaluate lots of different model options


# Predictions for 2017-2018 season - simple model, just Fourier terms, no covariates
# Train model based on prior data
train <- filter(flu_data_merge, season != "2017/2018" | week %in% 40:43) %>%
  # Create time series of ILI
  mutate(ILI = ts(ILI, frequency = 52, start = c(2006, 40)))

# Use auto.arima to determine best model
(fit <- auto.arima(train$ILI,
                   xreg = fourier(train$ILI, 3),
                   seasonal = FALSE,
                   lambda = -1))

# Fit ARIMA model based on most recent data using prior fit
pred_fit <- Arima(train$ILI, xreg = fourier(train$ILI, 3), model = fit)

# Simulate output
sim_output <- sample_predictive_trajectories_arima(
  pred_fit, 
  h = 28,
  xreg = fourier(train$ILI, 3, 28)
)

# Calculate forecast probabilities
forecast_results <- tibble() 

# Week ahead forecasts
for (i in 1:4) {
  for (j in seq(0, 12.9, 0.1)) {
    forecast_results <- bind_rows(
      forecast_results,
      tibble(target = paste(i, "wk ahead"),
             bin_start_incl = j,
             value = sum((sim_output[, i] == j))/nrow(sim_output)),
    )
  }
  forecast_results <- bind_rows(
    forecast_results,
    tibble(target = paste(i, "wk ahead"),
           bin_start_incl = 13,
           value = sum((sim_output[, i] >= 13))/nrow(sim_output))
  )
}

# Peak percentage forecasts
max_ili <- apply(sim_output, 1, max, na.rm = T)
for (j in seq(0, 12.9, 0.1)) {
  forecast_results <- bind_rows(
    forecast_results,
    tibble(target = "Season peak percentage",
           bin_start_incl = j,
           value = sum((max_ili == j))/length(max_ili))
  )
}
forecast_results <- bind_rows(
  forecast_results,
  tibble(target = "Season peak percentage",
         bin_start_incl = 13,
         value = sum((max_ili >= 13))/length(max_ili))
)

forecast_results <- forecast_results %>%
  mutate(location = "US National",
         type = "Bin",
         unit = "percent",
         bin_end_notincl = ifelse(bin_start_incl == 13, 100, 
                                  bin_start_incl + 0.1))

forecast(pred_fit,
         xreg=fourier(train$ILI, 3, h=31)) %>%
  autoplot()
?simulate.Arima
simulate(fit, nsim = 31, xreg = fourier(train$ILI, 3, h = 31))




snaive_gtrend <- snaive(flu_data_merge$ILI)
autoplot(snaive_gtrend)
checkresiduals(snaive_gtrend)
gglagplot(flu_data_merge$ILI)
ggsubseriesplot(flu_data_merge$ILI)
ggseasonplot(gtrend_US_ts)
?auto.arima
ggAcf(flu_data_merge$ILI)
diff(diff(flu_data_merge$ILI, 52), 1) %>% autoplot()
mstl(ili_US_ts) %>%
  forecast(method = "arima", bootstrap = TRUE) %>%
  autoplot()
ma2x52 <- ma(ili_US_ts, order = 52, centre = TRUE)
autoplot(ili_US_ts) +
  autolayer(ma2x52)
decompose(flu_data_merge$ILI) %>%
  autoplot()
seasonal::seas(ili_US_ts, x11 = "")
autoplot(gtrend_ILI_no2009$ILI)

BoxCox.lambda(ili_US_ts)

plots <- list()
for (i in 1:6) {
  fit <- auto.arima(flu_data_merge$ILI,
                      xreg = fourier(flu_data_merge$ILI, K = i),
                    seasonal = FALSE,
                    lambda = -1)
  plots[[paste(i)]] <- autoplot(forecast(fit,
                                         xreg = )) +
    xlab(paste("k=",i,"   AICC=",round(fit[["aicc"]],2))) +
    ylab("")
}

forecast(fit, xreg = fourier(gtrend_observed$ILI, K = 1, h = 30)) %>%
  autoplot()

checkresiduals(fit)
gridExtra::grid.arrange(
  plots[[1]],plots[[2]],plots[[3]],
  plots[[4]],plots[[5]],plots[[6]],
  nrow=3)

ggplot(data = gtrend_ILI_no2009, aes(x = hits, y = ILI)) +
  geom_point()
?simulate.Arima


fit <- auto.arima(gtrend_observed$ILI,
                  xreg = cbind(fourier(gtrend_observed$ILI, K = 4),
                               gtrend_observed$hits),
                  seasonal = FALSE,
                  lambda = -1)

test <- sample_predictive_trajectories_arima(
  fit,
  h = 52,
  xreg = cbind(fourier(gtrend_observed$ILI, K = 4, h = 52),
               snaive(ts(gtrend_observed$hits, frequency = 52, start = c(2006, 40)), h = 52)$mean)
)
?snaive()
hist(test[, 4])
dim(test)
(fit <- auto.arima(gtrend_observed$ILI, lambda = 0
                   xreg = cbind(fourier(ili_US_ts, 3),
                                gtrend_ILI_merge[1:618, c("hits", "h1_per_samples",
                                                       "h3_per_samples")]),
                   lambda = -1))
(fit)

checkresiduals(fit)

BoxCox(ili_US_ts, -1) %>% autoplot()

forecast(fit, 
         xreg = cbind(fourier(ili_US_ts, 3, 24),
                           snaive(ts(gtrend_ILI_merge$hits[1:723], frequency = 52, start = c(2004, 40)))$mean[1:24],
                           snaive(ts(gtrend_ILI_merge$h1_per_samples[1:723], frequency = 52, start = c(2004, 40)))$mean[1:24],
                           snaive(ts(gtrend_ILI_merge$h3_per_samples[1:723], frequency = 52, start = c(2004, 40)))$mean[1:24]),
         bootstrap = TRUE) %>%
  autoplot()

forecast(fit, level = 80) %>% autoplot()


## For future Google Trends - compare the 1 wk known value to what was observed at the same
# week last year. Adjust the 2, 3, 4 seasonal naive values by the ratio of the two 1 wk values