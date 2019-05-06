
Protea Analytics Influenza Forecasts
====================================

Prediction model for the CDC FluSight influenza forecasting challenge

This repo contains all documents related to the Protea Analytics entries for the US Centers for Disease Control and Prevention's 2018-2019 [FluSight influenza forecasting challenge](http://predict.cdc.gov), including conceptual documents, model fitting code, validation forecasts, and prospective forecasts submitted to CDC. Protea is also participating in the [FluSight Network](http://flusightnetwork.io/), a collaborative effort to build a weighted forecasting ensemble using data from multiple forecasting teams.

### Current forecasts

Forecasts are based on data from MMWR week 17, which encompasses Apr 21, 2019 to Apr 27, 2019. For interactive forecasts, please visit the [CDC](http://predict.cdc.gov) or [FluSight Network](http://flusightnetwork.io/) pages. The plots below illustrate the forecasts for weighted percentage of outpatient visits due to influenza-like illness (ILI) over the four weeks following the most recent data publication, with 80% prediction intervals.

<img src="README_files/figure-markdown_github/current forecasts-1.png" width="672" />

### Models

**Cheetah** - Cheetah is a weighted ensemble of the three other forecasts, based on leave-one-season-out cross validation of forecasts from the 2010-2011 through 2017-2018 seasons. Separate model weights are estimated for each month and target type (seasonal - onset, peak week, peak percentage; weekly - 1 to 4 weeks ahead). Weights were fit using a degenerate expectation maximization algorithm.

This model is being submitted to both CDC and the FluSight Network ensemble. Forecasts from this model can be found [here](CDC%20Submissions/2018-2019).

**Springbok** - Springbok is a dynamic harmonic regression model, using Fourier terms to capture the seasonality of the ILI trend along with ARIMA errors for short-term correlation. Region specific models include some combination of cumulative influenza subtype prevalence, national Google trends data, regional Google trends data, and estimates of ILI backfill as predictors in the model.

This model is being submitted to both CDC and the FluSight Network ensemble. Forecasts from this model can be found [here](CDC%20Submissions/2018-2019).

**Kudu** - Kudu is a historical average model weighted by the cumulative influenza A subtype prevelance. Separate Gaussian kernel density estimates based on previous ILI data were fit for H1N1 and H3N2 dominant seasons by weighting observed ILI values by the cumulative prevalence of that subtype in the season. Predictions are created by weighting those two kernel density estimates by the observed cumulative prevalence of the influenza A subtypes to date.

This model is being submitted to the FluSight Network ensemble only. Forecasts from this model can be found [here](Forecasts/2018-2019/Subtype%20Historical%20Average).

**Steenbok** - Steenbok is a simple historical average model. A Gaussian kernel density estimate based on previous ILI data is used to create predictions for the upcoming season.

This model is not submitted to either the CDC or FluSight Network and is only used as a component of the Cheetah model. Forecasts from this model can be found [here](Forecasts/2018-2019/Historical%20Average).
