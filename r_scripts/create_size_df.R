################################################################################
#' Author: Thomas Fjærvik
#' 
#' This script does the following:
#' 
#' Cleans the fund size data and saves the file as an Rdata object
################################################################################

################################################################################
#' Load packages
################################################################################

library(tidyverse)
library(readxl)

################################################################################
#' Load data
################################################################################

us_1_old <- read_excel("data/morningstar/aum_1_old.xlsx", skip = 8)
us_1_new <- read_excel("data/morningstar/aum_1.xlsx", skip = 8)
us_2_old <- read_excel("data/morningstar/aum_2_old.xlsx", skip = 8)
us_2_new <- read_excel("data/morningstar/aum_2.xlsx", skip = 8)
us_3_old <- read_excel("data/morningstar/aum_3_old.xlsx", skip = 8)
us_3_new <- read_excel("data/morningstar/aum_3.xlsx", skip = 8)
us_4_old <- read_excel("data/morningstar/aum_4_old.xlsx", skip = 8)
us_4_new <- read_excel("data/morningstar/aum_4.xlsx", skip = 8)


datasets <- list(us_1_old, us_1_new, us_2_old, us_2_new, us_3_old, us_3_new, us_4_old, us_4_new)

# delete first two rows of all datasets

datasets <- lapply(datasets, function(x) x[-c(1:2), ])

# change columns names
# # Put data on long format
# deselect `Group/Investment`
# fix date format
# create YearMonth variable
# deselect date

datasets <- map(datasets, function(x) rename_with(x, 
                ~sub("Fund Size - comprehensive \\(Monthly\\) ", "", colnames(x))
                )) %>% 
  map(., ~pivot_longer(.x, 
                       cols = !c(1:3),
                       names_to = "date",
                       values_to = "size")) %>% 
  map(., ~{.x %>% 
      select(FundId, SecId, date, size) %>% 
      mutate(date = sub("\\.", "-", date),
      date = paste0("01-", date),
      date = base::as.Date(date, "%d-%m-%Y"),
      YearMonth = zoo::as.yearmon(date)
      ) %>% 
      select(-date)
    })


# Create one dataset and arrange by FundId, SecId and Yearmonth

size <- datasets %>% 
  data.table::rbindlist() %>% 
  arrange(FundId, SecId, YearMonth)

# Save file

save(size, file = "pipeline/create_size_df/size.Rdata")


