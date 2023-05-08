################################################################################
#' Author: Thomas Fjærvik
#' 
#' This script does the following:
#' 
#' Merges ltd estimates with Fama-French data, the momentum factor and market value
#' data and creates return-based variables for subsequent analysis.
################################################################################

################################################################################
#' Load packages
################################################################################

library(tidyverse)
library(readxl)
library(readr)

################################################################################
#' Load data
################################################################################

# LTD estimates
load("pipeline/read_rds_data/ltd_estimates_monthly.Rdata")

# Fama-French factors and risk-free rate
ff <- read_csv("data/kf_data_library/ff3_monthly.csv")
colnames(ff) <- c("date", "excess_ret", "smb", "hml", "rf")

ff <- ff %>% 
  mutate(date = paste0(substr(date, 1, 4), "-", substr(date, 5, 6), "-01"),
         date = as.Date(date), 
         YearMonth = zoo::as.yearmon(date),
         excess_ret = excess_ret/100,
         smb = smb/100,
         hml = hml/100,
         rf = rf/100) %>% 
  select(YearMonth, smb, hml, rf)

# Momentum factor 

mom <- read_csv("data/kf_data_library/momentum_monthly.csv")

mom <- mom %>% 
  rename(date = Date) %>% 
  mutate(date = paste0(substr(date, 1, 4), "-", substr(date, 5, 6), "-01"),
         date = as.Date(date), 
         YearMonth = zoo::as.yearmon(date),
         mom = mom/100) %>% 
  select(YearMonth, mom)

# Merge the two factor datasets

factors <- ff %>% 
  left_join(mom, by = "YearMonth")

# Market value

load("pipeline/create_size_df/size.Rdata")

################################################################################
#' Merge ltd estimates with factors and market value data
################################################################################

# Merge with factors 

ltd_estimates_monthly <- ltd_estimates_monthly %>% 
  left_join(factors, by = "YearMonth")

# Calculate excess returns

ltd_estimates_monthly <- ltd_estimates_monthly %>% 
  mutate(xret = monthly_ret - rf,
         xret_index = monthly_ret_index - rf) %>% 
  select(-c(ret, ret_index, monthly_estimate))

# Merge with market value data

ltd_estimates_monthly <- ltd_estimates_monthly %>% 
  left_join(size, by = c("FundId", "SecId", "YearMonth"))

################################################################################
#' Save dataset
################################################################################

save(ltd_estimates_monthly, file = "pipeline/merge_ltd_estimates_with_factors/monthly_df_before_return_based_variables.Rdata")
