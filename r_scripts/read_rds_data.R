################################################################################
#' Author: Thomas Fjærvik
#' 
#' This script does the following:
#' 
#' Read processed RDS-files into R and create a dataset ready for 
#' calculateing return based variables
################################################################################

################################################################################
#' Load packages
################################################################################


library(tidyverse)
library(zoo)
library(furrr)
library(readxl)
library(lubridate)
library(data.table)

################################################################################
#' Load RDS data
################################################################################

ltd_estimates <-  list.files(path = "data/ltd_estimates", 
                             pattern = ".rds", full.names = TRUE) %>%
  map(readRDS) %>% 
  #map(~{.x %>% 
  #      drop_na(monthly_estimate)}
  #   ) %>% 
  rbindlist()

# Look for duplicated dates

ltd_estimates %>% 
  group_by(FundId) %>% 
  summarise(unique = length(unique(date)), total = length(date)) %>% 
  mutate(diff = total-unique,
         index = diff != 0) %>% 
  pull(index) %>% 
  which()

# There are no duplicated dates

################################################################################
#' Split the estimates into different columns
################################################################################

ltd_estimates <- ltd_estimates %>% 
  mutate(LTD_estimate = str_replace_all(monthly_estimate, "(?<!e)-", "/")) %>% 
  separate("LTD_estimate", c("ltd", "cop_combo", "utd", "ltd_nll", "cop_combo_nll"),
           sep = "/")

# Make the new columns numeric

ltd_estimates <- ltd_estimates %>% 
  mutate(across(c(ltd, cop_combo, utd, ltd_nll, cop_combo_nll), as.numeric))

# Need to add NAV to the funds so we can calculate value-weighted returns.
# In the meantime, we save the dataset before removing NA values, so that we may
# calculate value-weighted returns. First, we aggregate returns to a monthly level.
# Must first divide by 100

ltd_estimates_monthly <- ltd_estimates %>% 
  mutate(ret_index = as.numeric(ret_index),
         ret = ret/100,
         ret_index = ret_index/100) %>% 
  group_by(FundId, year, month) %>% 
  mutate(monthly_ret = prod(1+ret) - 1,
         monthly_ret_index = prod(ret_index + 1) - 1) %>% 
  filter(date == max(date)) %>% 
  ungroup()


# Save monthly dataset

save(ltd_estimates_monthly, file = "pipeline/read_rds_data/ltd_estimates_monthly.Rdata")

