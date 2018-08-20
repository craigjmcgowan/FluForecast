# Functions to keep weeks in order and reset them
week_inorder <- function(week, season){
  case_when(as.numeric(week) < 40 & 
              season %in% c("1997/1998", "2003/2004", "2008/2009", "2014/2015") ~ 
              as.numeric(week) + 53,
            as.numeric(week) < 40 ~ as.numeric(week) + 52,
            TRUE ~ as.numeric(week))
}

week_reset <- function(week, season) {
  case_when(as.numeric(week) > 53 & 
              season %in% c("1997/1998", "2003/2004", "2008/2009", "2014/2015") ~ 
              as.numeric(week) - 53,
            as.numeric(week) > 52 & 
              !season %in% c("1997/1998", "2003/2004", "2008/2009", "2014/2015") ~
              as.numeric(week) - 52,
            TRUE ~ as.numeric(week))
}
