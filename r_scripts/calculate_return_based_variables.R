################################################################################
#' Author: Thomas Fjærvik
#' 
#' This script does the following:
#' 
#' Calculates return based variables to be used in subsequent analysis
################################################################################

################################################################################
#' Load packages
################################################################################

library(tidyverse)
library(data.table)


################################################################################
#' Load data
################################################################################

load("pipeline/merge_ltd_estimates_with_factors/monthly_df_before_return_based_variables.Rdata")

################################################################################
#' Load functions
################################################################################

source("r_scripts/function_scripts/return_based_variables_functions.R")

################################################################################
#' Calculate return based variables
################################################################################

# Split into list and operate on each dataset separately. data.table::rbindlist()
# at the end. 

data_list <- split(ltd_estimates_monthly, f = ltd_estimates_monthly$FundId)

# Make all the elements tibbles

data_list <- data_list %>% 
  map(., as_tibble)


# test_data_list <- data_list %>% 
#   map(., ~{.x %>% 
#       slice_head(n = 15)}) %>% 
#   head(n = 10)

return_based_variables <- data_list %>% 
  map(., ~ {rolling_window_multivariate_modified(data = .x,
                                                 columns = c("xret", "xret_index"),
                                                 function_name = calc_beta,
                                                 output_variable = "beta")}) %>% 
  map(., ~ {rolling_window_multivariate_modified(data = .x,
                                                 columns = c("xret", "xret_index"),
                                                 function_name = calc_cokurtosis,
                                                 output_variable = "cokurtosis")}) %>% 
  map(., ~ {rolling_window_multivariate_modified(data = .x,
                                                 columns = c("xret", "xret_index"),
                                                 function_name = calc_coskewness,
                                                 output_variable = "coskewness")}) %>% 
  map(., ~ {rolling_window_multivariate_modified(data = .x,
                                                 columns = c("xret", "xret_index"),
                                                 function_name = calc_downside_beta,
                                                 output_variable = "downside_beta")}) %>% 
  map(., ~ {rolling_window_multivariate_modified(data = .x,
                                                 columns = c("xret", "xret_index"),
                                                 function_name = calc_upside_beta,
                                                 output_variable = "upside_beta")}) %>% 
  map(., ~ {rolling_window_multivariate_modified(data = .x,
                                                 columns = c("xret", "xret_index"),
                                                 function_name = capm_alpha,
                                                 output_variable = "capm_alpha")}) %>% 
  map(., ~ {rolling_window_multivariate_modified(data = .x,
                                                 columns = c("xret", "xret_index", "smb", "hml", "mom"),
                                                 function_name = car_alpha,
                                                 output_variable = "car_alpha")}) %>% 
  map(., ~ {rolling_window_multivariate_modified(data = .x,
                                                 columns = c("xret", "xret_index", "smb", "hml"),
                                                 function_name = ff3_alpha,
                                                 output_variable = "ff3_alpha")})

# Create one dataset

fund_data <- return_based_variables %>% 
  data.table::rbindlist() %>% 
  as_tibble() %>% 
  arrange(FundId, SecId, YearMonth)

# Save dataset

save(fund_data, file = "pipeline/calculate_return_based_variables/fund_data.Rdata")


