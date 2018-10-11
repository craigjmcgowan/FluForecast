# Create ensemble weights

# This code borrows heavily from the estimation code used in the FluSight Network
#   efforts to create a weighted forecasting ensemble. Code can be found on GitHub
#   at https://github.com/FluSightNetwork/cdc-flusight-ensemble

library(tidyverse)
library(FluSight)
library(pipeR)
library(parallel)

# Load scores
load("Data/model_scores.Rdata")
source("R/utils.R")
source("R/degenerate_em_functions.R")

# Set up cores
options(mc.cores=parallel::detectCores()-1L)

# Define weight structures
#   Equal weights
#   Constant weights
#   Week weights
#   4-week weights
#   Region weights
#   Target-type weights
#   Target weights
#   Week, region weights
#   Week, target-type weights
#   Week, target weights
#   4-week, region weights
#   4-week, target-type weights
#   4 week, target weights
#   Region, target-type weights
#   Region, target weights
#   Week, region, target-type  weights
#   Week, region, target weights
#   4-week, region, target-type  weights
#   4-week, region, target weights

# Create equal weights ------
nteams <- length(unique(all_eval_scores$team[!grepl('ens', all_eval_scores$team)]))

equal_weights_df <- crossing(component_model_id = 
                               unique(all_eval_scores$team[!grepl('ens', all_eval_scores$team)]),
                             season = c(unique(all_eval_scores$season),
                                        paste0(substr(last(unique(all_eval_scores$season)),
                                                      6, 9), "/",
                                               as.numeric(substr(last(unique(all_eval_scores$season)),
                                                      6, 9)) + 1)),
                             target = unique(all_eval_scores$target)) %>%
  mutate(weight = 1 / nteams)

write.csv(equal_weights_df, "weights/equal-weights.csv", row.names = FALSE, quote=FALSE)

## Specify target types:
target.types.list = list(
  "seasonal"=c("Season onset", "Season peak week", "Season peak percentage"),
  "weekly"=paste0(1:4," wk ahead")
)
## Specify 'months'
month.types.list = list(
  "October" = as.character(40:44),
  "November" = as.character(45:48),
  "December" = as.character(49:52),
  "January" = as.character(53:56),
  "February" = as.character(57:60),
  "March" = as.character(61:64),
  "April" = as.character(65:68),
  "May" = as.character(69:72)
)

## Specify portions of cv_apply indexer lists corresponding to 
## Model Week, Location, Target:
weighting.scheme.partial.indexer.lists = list(
  "constant" = list(all = NULL, all = NULL, all = NULL),
  "region-based" = list(all = NULL, each = NULL, all = NULL),
  "target-type-based" = list(all = NULL, all = NULL, 
                             subsets = target.types.list),
  "target-based" = list(all = NULL, all = NULL, each = NULL),
  "region-target-type-based" = list(all = NULL, each = NULL, 
                                    subsets = target.types.list),
  "region-target-based" = list(all = NULL, each = NULL, each = NULL)
)

week.partial.indexer.lists <- list(
  "week-based" = list(each = NULL, all = NULL, all = NULL),
  "month-based" = list(subsets = month.types.list, all = NULL, all = NULL),
  "week-region-based" = list(each = NULL, each = NULL, all = NULL),
  "week-target-type-based" = list(each = NULL, all = NULL, 
                                  subsets = target.types.list),
  "week-target-based" = list(each = NULL, all = NULL, each = NULL),
  "month-region-based" = list(subsets = month.types.list, each = NULL, all = NULL),
  "month-target-type-based" = list(subsets = month.types.list, all = NULL, 
                                  subsets = target.types.list),
  "month-target-based" = list(subsets = month.types.list, all = NULL, each = NULL),
  "week-target-type-region-based" = list(each = NULL, each = NULL,
                                         subsets = target.types.list),
  "week-target-region-based" = list(each = NULL, each = NULL, each = NULL),
  "month-target-type-region-based" = list(subsets = month.types.list, each = NULL,
                                         subsets = target.types.list),
  "month-target-region-based" = list(subsets = month.types.list,
                                      each = NULL, each = NULL)
)

# Set up variable with scores grouped by necessary variables
cv_full_scores <- all_full_scores %>%
  # Remove ensemble models
  filter(!grepl('ens', team), !grepl('Ens', team)) %>%
  # Create ordered week variable
  mutate(order_week = case_when(season == "2014/2015" & forecast_week < 40 ~ 
                                  forecast_week + 53,
                                forecast_week < 40 ~ forecast_week + 52,
                                TRUE ~ forecast_week),
         metric = "some log score",
         score_to_optimize = score) %>%
  # Remove above week 71 - never in play for scoring
  filter(order_week != 73) %>%
  # Cast to array for weight estimating
  reshape2::acast(season ~ order_week ~ location ~ target ~ metric ~ team, 
                  value.var="score_to_optimize") %>>%
                  {names(dimnames(.)) <- c("Season", "Model Week", "Location", "Target", "Metric", "Model"); .}

cv_eval_scores <- all_eval_scores %>%
  # Remove ensemble models
  filter(!grepl('ens', team), !grepl('Ens', team)) %>%
  # Create ordered week variable
  mutate(order_week = case_when(season == "2014/2015" & forecast_week < 40 ~ 
                                  forecast_week + 53,
                                forecast_week < 40 ~ forecast_week + 52,
                                TRUE ~ forecast_week),
         metric = "some log score",
         score_to_optimize = score) %>%
  # Remove week 71 - never in play for scoring
  filter(order_week != 71) %>%
  # Cast to array for weight estimating
  reshape2::acast(season ~ order_week ~ location ~ target ~ metric ~ team, 
                  value.var="score_to_optimize") %>>%
  {names(dimnames(.)) <- c("Season", "Model Week", "Location", "Target", "Metric", "Model"); .}
  
## Indexer lists for prospective forecasts:
weighting.scheme.prospective.indexer.lists =
  weighting.scheme.partial.indexer.lists %>>%
  lapply(function(partial.indexer.list) {
    c(list(all=NULL), # Season, Model Week
      partial.indexer.list, # Location, Target
      list(all=NULL), # Metric
      list(all=NULL) # Model should always be all=NULL
    )
  })

## Indexer lists for CV forecasts:
weighting.scheme.cv.indexer.lists =
  weighting.scheme.partial.indexer.lists %>>%
  lapply(function(partial.indexer.list) {
    c(list(loo=NULL), # Season, Model Week
      partial.indexer.list, # Location, Target
      list(all=NULL), # Metric
      list(all=NULL) # Model should always be all=NULL
    )
  })


## Indexer lists for weekly weighted prospective forecasts:
week.prospective.indexer.lists =
  week.partial.indexer.lists %>>%
  lapply(function(partial.indexer.list) {
    c(list(all=NULL), # Season, Model Week
      partial.indexer.list, # Location, Target
      list(all=NULL), # Metric
      list(all=NULL) # Model should always be all=NULL
    )
  })

## Indexer lists for weekly CV forecasts:
week.cv.indexer.lists =
  week.partial.indexer.lists %>>%
  lapply(function(partial.indexer.list) {
    c(list(loo=NULL), # Season, Model Week
      partial.indexer.list, # Location, Target
      list(all=NULL), # Metric
      list(all=NULL) # Model should always be all=NULL
    )
  })

## Target-type df:
target.types.df =
  lapply(target.types.list, tibble::as_tibble) %>>%
  dplyr::bind_rows(.id="target.type") %>>%
  dplyr::rename(target=value)

month.types.df = 
  lapply(month.types.list, tibble::as_tibble) %>>%
  dplyr::bind_rows(.id="month") %>>%
  dplyr::rename(`Model Week`=value)


#### Generate the weight files: #####

## Non-weekly weights
for (weighting.scheme.i in seq_along(weighting.scheme.partial.indexer.lists)) {
  ## extract info from lists:
  weighting.scheme.name = names(weighting.scheme.partial.indexer.lists)[[weighting.scheme.i]]
  weighting.scheme.cv.indexer.list = weighting.scheme.cv.indexer.lists[[weighting.scheme.i]]
  weighting.scheme.prospective.indexer.list = weighting.scheme.prospective.indexer.lists[[weighting.scheme.i]]
  print(weighting.scheme.name)
  ## determine season label for next season:
  cv.season.ints = as.integer(gsub("/.*","",dimnames(cv_eval_scores)[["Season"]]))
  prospective.season.int = max(cv.season.ints) + 1L
  prospective.season.label = prospective.season.int %>>% paste0("/",.+1L)
  ## generate weight df's and bind together:
  cv.weight.df = generate_indexer_list_weights(
    cv_eval_scores, weighting.scheme.cv.indexer.list
  )
  prospective.weight.df = generate_indexer_list_weights(
    cv_eval_scores, weighting.scheme.prospective.indexer.list
  ) %>>%
    dplyr::mutate(season=prospective.season.label)
  combined.weight.df =
    dplyr::bind_rows(cv.weight.df, prospective.weight.df)
  ## expand out target types if applicable:
  if ("target" %in% names(combined.weight.df)) {
    combined.weight.df <-
      combined.weight.df %>>%
      dplyr::rename(target.type=target) %>>%
      dplyr::left_join(
        dimnames(cv_eval_scores)[["Target"]] %>>%
        {tibble::tibble(target.type=., target=.)} %>>%
          dplyr::bind_rows(target.types.df),
        by = "target.type"
      ) %>>%
      dplyr::select(-target.type) %>>%
      ## restore original column order:
      magrittr::extract(names(combined.weight.df))
  }
  ## print weight df and write to csv file:
  print(combined.weight.df)
  write.csv(combined.weight.df, file.path("weights",paste0(weighting.scheme.name,"-weights.csv")), row.names = FALSE, quote=FALSE)
}


## Weekly weights
for (weighting.scheme.i in seq_along(week.partial.indexer.lists)) {
  ## extract info from lists:
  weighting.scheme.name = names(week.partial.indexer.lists)[[weighting.scheme.i]]
  weighting.scheme.cv.indexer.list = week.cv.indexer.lists[[weighting.scheme.i]]
  weighting.scheme.prospective.indexer.list = week.prospective.indexer.lists[[weighting.scheme.i]]
  print(weighting.scheme.name)
  ## determine season label for next season:
  cv.season.ints = as.integer(gsub("/.*","",dimnames(cv_eval_scores)[["Season"]]))
  prospective.season.int = max(cv.season.ints) + 1L
  prospective.season.label = prospective.season.int %>>% paste0("/",.+1L)
  ## generate weight df's and bind together:
  cv.weight.df = generate_indexer_list_weights(
    cv_full_scores, weighting.scheme.cv.indexer.list
  )
  prospective.weight.df = generate_indexer_list_weights(
    cv_full_scores, weighting.scheme.prospective.indexer.list
  ) %>>%
    dplyr::mutate(season=prospective.season.label)
  combined.weight.df =
    dplyr::bind_rows(cv.weight.df, prospective.weight.df)
  ## expand out target types if applicable:
  if ("target" %in% names(combined.weight.df)) {
    combined.weight.df <-
      combined.weight.df %>>%
      dplyr::rename(target.type=target) %>>%
      dplyr::left_join(
        dimnames(cv_full_scores)[["Target"]] %>>%
        {tibble::tibble(target.type=., target=.)} %>>%
          dplyr::bind_rows(target.types.df),
        by = "target.type"
      ) %>>%
      dplyr::select(-target.type) %>>%
      ## restore original column order:
      magrittr::extract(names(combined.weight.df))
  }
  ## expand out months into weeks if applicable:
  if ("Model Week" %in% names(combined.weight.df)) {
    combined.weight.df <-
      combined.weight.df %>>%
      dplyr::rename(month=`Model Week`) %>>%
      dplyr::mutate(month = as.character(month)) %>>%
      dplyr::left_join(
        dimnames(cv_full_scores)[["Model Week"]] %>>%
        {tibble::tibble(month=., `Model Week`=.)} %>>%
          dplyr::bind_rows(month.types.df),
        by = "month"
      ) %>>%
      dplyr::select(-month) %>>%
      ## restore original column order:
      magrittr::extract(names(combined.weight.df))
  }
  ## print weight df and write to csv file:
  print(combined.weight.df)
  write.csv(combined.weight.df, file.path("weights",paste0(weighting.scheme.name,"-weights.csv")), row.names = FALSE, quote=FALSE)
}

