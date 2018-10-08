# gtrendsR functions modified to ignore cookies
require(gtrendsR)

get_api_cookies_hack <- function (cookie_url) 
{
  cookie_handler <- curl::new_handle()
  cookie_req <- curl::curl_fetch_memory(cookie_url, handle = cookie_handler)
  curl::handle_cookies(cookie_handler)
  return(cookie_handler)
}

get_widget_hack <- function (comparison_item, category, gprop, hl, cookie_url, 
          tz, cookie_handler) 
{
  token_payload <- list()
  token_payload$comparisonItem <- comparison_item
  token_payload$category <- category
  token_payload$property <- gprop
  url <- URLencode(paste0("https://www.google.com/trends/api/explore?property=&req=", 
                          jsonlite::toJSON(token_payload, auto_unbox = TRUE), 
                          "&tz=", tz, "&hl=", hl))
  url <- gtrendsR:::encode_keyword(url)
  if (is.null(cookie_handler)) {
    cookie_handler <- get_api_cookies_hack(cookie_url)
  }
  widget <- curl::curl_fetch_memory(url, handle = cookie_handler)
  stopifnot(widget$status_code == 200)
  temp <- rawToChar(widget$content)
  Encoding(temp) <- "UTF-8"
  myjs <- jsonlite::fromJSON(substring(temp, first = 6))
  widget <- myjs$widgets
}

interest_over_time_hack <- function (widget, comparison_item, tz, cookie_handler) 
{
  payload2 <- list()
  if (!is.na(widget$request$locale[1])) {
    payload2$locale <- widget$request$locale[1]
    payload2$comparisonItem <- widget$request$comparisonItem[[1]]
    payload2$resolution <- widget$request$resolution[1]
    payload2$requestOptions$category <- widget$request$requestOptions$category[1]
    payload2$requestOptions$backend <- widget$request$requestOptions$backend[1]
    payload2$time <- widget$request$time[1]
    payload2$requestOptions$property <- widget$request$requestOptions$property[1]
    token_payload2 <- widget$token[1]
  }
  else {
    payload2$locale <- widget$request$locale[2]
    payload2$comparisonItem <- widget$request$comparisonItem[[2]]
    payload2$resolution <- widget$request$resolution[2]
    payload2$requestOptions$category <- widget$request$requestOptions$category[2]
    payload2$requestOptions$backend <- widget$request$requestOptions$backend[2]
    payload2$time <- widget$request$time[2]
    payload2$requestOptions$property <- widget$request$requestOptions$property[2]
    token_payload2 <- widget$token[2]
  }
  url <- URLencode(paste0("https://www.google.com/trends/api/widgetdata/multiline/csv?req=", 
                          jsonlite::toJSON(payload2, auto_unbox = T), "&token=", 
                          token_payload2, "&tz=", tz))
  url <- gtrendsR:::encode_keyword(url)
  res <- curl::curl_fetch_memory(url, handle = cookie_handler)
  stopifnot(res$status_code == 200)
  con <- textConnection(rawToChar(res$content))
  df <- read.csv(con, skip = 1, stringsAsFactors = FALSE)
  close(con)
  if (nrow(df) < 1) {
    return(NULL)
  }
  n <- nrow(df)
  df <- reshape(df, varying = names(df)[2:ncol(df)], v.names = "hits", 
                direction = "long", timevar = "temp", times = names(df)[2:ncol(df)])
  df$temp <- NULL
  df <- cbind(df, comparison_item[rep(seq_len(nrow(comparison_item)), 
                                      each = n), 1:2], row.names = NULL)
  df$geo <- ifelse(df$geo == "", "world", df$geo)
  df$gprop <- ifelse(widget$request$requestOptions$property[1] == 
                       "", "web", widget$request$requestOptions$property[1])
  df$category <- widget$request$requestOptions$category[1]
  names(df)[1] <- "date"
  df$id <- NULL
  if (unique(comparison_item$time) == "all") {
    df$date <- anytime::anydate(df$date)
  }
  else {
    df$date <- anytime::anytime(df$date)
  }
  return(df)
}


gtrends_hack <- function (keyword = NA, geo = "", time = "today+5-y", gprop = c("web", 
                                                                "news", "images", "froogle", "youtube"), category = 0, hl = "en-US", 
          low_search_volume = FALSE, cookie_url = "http://trends.google.com/Cookies/NID", 
          TZ = -120, cookie_handler = NULL) 
{
  stopifnot((length(keyword)%%length(geo) == 0) || (length(geo)%%length(keyword) == 
                                                      0), is.vector(keyword), length(keyword) <= 5, length(geo) <= 
              5, length(time) == 1, length(hl) == 1, is.character(hl), 
            length(cookie_url) == 1, 
            is.character(cookie_url))
  
  gprop <- match.arg(gprop, several.ok = FALSE)
  gprop <- ifelse(gprop == "web", "", gprop)
  comparison_item <- data.frame(keyword, geo, time, stringsAsFactors = FALSE)
  widget <- get_widget_hack(comparison_item, category, gprop, hl, 
                       cookie_url, TZ, cookie_handler)
  interest_over_time2 <- interest_over_time_hack(widget, comparison_item, 
                                           TZ, cookie_handler)
  interest_by_region <- gtrendsR:::interest_by_region(widget, comparison_item, 
                                           low_search_volume, TZ)
  related_topics <- gtrendsR:::related_topics(widget, comparison_item, 
                                   hl, TZ)
  related_queries <- gtrendsR:::related_queries(widget, comparison_item, 
                                     TZ)
  res <- list(interest_over_time = interest_over_time, interest_by_country = do.call(rbind, 
                                                                                     interest_by_region[names(interest_by_region) == "country"]), 
              interest_by_region = do.call(rbind, interest_by_region[names(interest_by_region) == 
                                                                       "region"]), interest_by_dma = do.call(rbind, interest_by_region[names(interest_by_region) == 
                                                                                                                                         "dma"]), interest_by_city = do.call(rbind, interest_by_region[names(interest_by_region) == 
                                                                                                                                                                                                         "city"]), related_topics = related_topics, related_queries = related_queries)
  res <- lapply(res, function(x) {
    row.names(x) <- NULL
    x
  })
  class(res) <- c("gtrends", "list")
  return(res)
}
rm(this_handle, cookie_handler)
cookie_handler <- this_handle <- get_api_cookies_hack("http://trends.google.com/Cookies/NID")
MA_flu_3 <- gtrends_hack(keyword = "influenza",
                         geo = "US-MA",
                         cookie_handler = this_handle)$interest_over_time

MA_flu_1 <- gtrends(keyword = "flu", geo = "US-MA")$interest_over_time
