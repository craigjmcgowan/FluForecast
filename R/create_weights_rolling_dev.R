# Create ensemble weights

# This code borrows heavily from the estimation code used in the FluSight Network
#   efforts to create a weighted forecasting ensemble. Code can be found on GitHub
#   at https://github.com/FluSightNetwork/cdc-flusight-ensemble

library(tidyverse)
library(FluSight)
library(pipeR)
library(parallel)

source("R/utils.R")
source("R/degenerate_em_functions.R")

# Load scores - use for completed seasons
# weekly_scores_1819 <- readRDS('Data/weekly_scores_1819.Rds')
# weekly_eval_scores_1819 <- readRDS('Data/weekly_eval_scores_1819.Rds')
# weekly_scores_1920 <- readRDS('Data/weekly_scores_1920.Rds')
# weekly_eval_scores_1920 <- readRDS('Data/weekly_eval_scores_1920.Rds')

# Set most recent week that data is available
this_week <- 'EW44'

# Load relevant scores from the current season
all_full_scores <- bind_rows(readRDS('Data/weekly_scores_1920.Rds')[[this_week]])


# Set up cores
options(mc.cores=parallel::detectCores()-1L)

# Define prospective season
pro_season <- "2019/2020"

# Define weight structures
#   Constant weights
#   Region weights
#   Target-type weights
#   Target weights
#   Region, target-type weights
#   Region, target weights

# Create equal weights ------
nteams <- length(unique(all_eval_scores$team[!grepl('ens', all_eval_scores$team)]))
cv_seasons <- unique(all_eval_scores$season)[-which(unique(all_eval_scores$season) == pro_season)]

# equal_weights_df <- crossing(component_model_id = 
#                                unique(all_eval_scores$team[!grepl('ens', all_eval_scores$team)]),
#                              season = c(cv_seasons, pro_season),
#                              target = unique(all_eval_scores$target)) %>%
#   mutate(weight = 1 / nteams)
# 
# write.csv(equal_weights_df, "weights/equal-weights.csv", row.names = FALSE, quote=FALSE)

## Specify target types:
target.types.list = list(
  "seasonal" = c("Season onset", "Season peak week", "Season peak percentage"),
  "weekly" = paste0(1:4," wk ahead")
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

## Indexer lists for ongoing forecasts:
weighting.scheme.ongoing.indexer.lists =
  weighting.scheme.partial.indexer.lists %>>%
  lapply(function(partial.indexer.list) {
    c(list(all=NULL), # Season - there should only be one season in these data 
      partial.indexer.list, # Model week, Location, Target
      list(all=NULL), # Metric - only one metric, so just leave as all = NULL
      list(all=NULL) # Model should always be all=NULL - include all models for potential weights
    )
  })

# Set up variable with scores grouped by necessary variables
cv_full_scores <- all_full_scores %>%
  # Remove prospective season scores if they exist
  filter(season == pro_season) %>%
  # Remove ensemble models
  filter(!grepl('ens', team), !grepl('Ens', team)) %>%
  # Create ordered week variable
  mutate(order_week = week_inorder(forecast_week, season),
         metric = "some log score",
         score_to_optimize = score) %>%
  # Remove order week 73 from seasons that have 53 MMWR weeks - never in play for scoring
  filter(order_week != 73) %>%
  # Cast to array for weight estimating
  reshape2::acast(season ~ order_week ~ location ~ target ~ metric ~ team, 
                  value.var="score_to_optimize") %>>%
  {names(dimnames(.)) <- c("Season", "Model Week", "Location", "Target", "Metric", "Model"); .}


## Target-type df:
target.types.df =
  lapply(target.types.list, tibble::as_tibble) %>>%
  dplyr::bind_rows(.id="target.type") %>>%
  dplyr::rename(target=value)


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
  prospective.season.int = max(cv.season.ints) #+ 1L
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
  write.csv(combined.weight.df, 
            file.path("weights",
                      paste0(weighting.scheme.name,"-weights.csv")),
            row.names = FALSE, quote=FALSE)
}



