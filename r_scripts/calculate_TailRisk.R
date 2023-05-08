################################################################################
#' Author: Thomas Fjærvik
#' 
#' This script does the following:
#' 
#' Calculates the TailRisk measure from the paper "Tail risk in hedge funds:
#' A unique view from portfolio holdings" (Agarwal, Ruenzi, Weigert 2017). The 
#' measure consists of the two measures TailSens and Expected Shortfall
################################################################################

################################################################################
#' Load packages
################################################################################

library(tidyverse)

################################################################################
#' Load data
################################################################################

load("pipeline/calculate_return_based_variables/fund_data.Rdata")

################################################################################
#' Calculate TailSens
################################################################################

# We use rolling windows of 24 months. The estimation is nonparametric and based
# on the empirical return distribution of fund i and the corresponding market 
# return. We use a cutoff of q=0.05 (quantile).

test_vec <- seq(1:15)
test_vec_2 <- fund_data$YearMonth %>% head(n=15)
# test_vec_2 <- seq(101, 115, by = 1)

test_tibble <- tibble(
  data = test_vec,
  index = test_vec_2)


test_tibble <- test_tibble %>% 
  mutate(rolling_mean = slider::slide_index(data, .i = index, .f = ~mean(.x), 
                                            .before = 1, 
                                            .after = -1/12,
                                            .complete = TRUE))

test_tibble_2 <- test_tibble %>% 
  mutate(rolling_mean = slider::slide_dbl(data, mean, 
                                          .before = 12, 
                                          .after = -1,
                                          .complete = TRUE))


# We can use this method to compute TailSens after we have written a function for 
# calculating TailSens

# x <- test_data
# market_column <- "xret_index"
# ret_colum <- "xret"

TailSens <- function(x, market_column, ret_colum){
  
  # x is the dataset
  # columns is the columns for which to calculate TailSens. Must be written
  # as character names
  
  cutoff_market <- quantile(x[[market_column]], probs = 0.05, type = 1)
  cutoff_ret <- quantile(x[[ret_colum]], probs = 0.05, type = 1)
  
  market_observations <- x %>% 
    filter(xret_index <= cutoff_market)
  
  length_market <- nrow(market_observations)
  
  market_intersect_ret_observations <- market_observations %>% 
    filter(xret <= cutoff_ret)
  
  length_ret <- nrow(market_intersect_ret_observations)
  
  ltd <- length_ret/length_market
  
  # ltd <- x %>% 
  #   filter(xret_index <= cutoff) %>% 
  #   select(xret) %>% 
  #   summarise(ltd = mean(xret)) %>% 
  #   as.numeric()
  
  return(ltd)
  
} 

# test_data <- fund_data %>% head(n=24)
# test_data %>% TailSens(., market_column = "xret_index", ret_colum = "xret")

# test if it is possible to get 0.5 TailSens

test_vec <- c(rep(1, times= 22), 0, -1)
test_vec_2 <- c(-1, rep(1, times = 22), 0)
test_vec_3 <- fund_data$YearMonth %>% head(n=24)

test_data <- tibble(
  YearMonth = test_vec_3,
  xret = test_vec_2,
  xret_index = test_vec
)

TailSens(test_data, "xret_index", "xret") # Yes, it is possible
# Try to compute TailSens for the whole dataset

fund_data_2 <- fund_data %>% 
  #head(n = 10000) %>% 
  group_by(FundId) %>% 
  nest() %>% 
  mutate(TailSens_list = map(data, ~ .x %>% 
                               mutate(TailSens = slider::slide_index(.x , .i = YearMonth, .f = ~TailSens(., 
                                                                          "xret_index", 
                                                                          "xret"), 
                             .before = 2, # 24 months
                             .after = -1/12, # 1 month
                             .complete = TRUE))
  )
  ) %>% 
  select(-data) %>% 
  unnest(TailSens_list)

# Save file

save(fund_data_2, file = "pipeline/calculate_TailRisk/TailSens.Rdata")

################################################################################
#' Compute Expected shortfall
################################################################################

library(PerformanceAnalytics)

test_data <- fund_data %>% head(n=24)

es <- ES(test_data$xret, p = 0.95, method = "historical") %>% as.numeric()

# Now test by doing it manually.

test_data_es <- test_data %>%
  mutate(cutoff = quantile(xret, probs = 0.05)) %>%
  filter(xret <= cutoff) %>%
  summarise(es = mean(xret))

# Write function to calculate expected shortfall to be used in rolling windows

expected_shortfall <- function(x, column){
  
  # x is the data
  # column is the column for which to calculate expected shortfall. Must be
  # a character
  
  PerformanceAnalytics::ES(x[[column]], p = 0.95, method = "historical") %>% as.numeric()
  
}

# Calculate expected shortfall both for market return and fund returns by group

fund_data_3 <- fund_data_2 %>% 
  group_by(FundId) %>% 
  nest() %>% 
  mutate(ES_list = map(data, ~ .x %>% 
                               mutate(ES_market = slider::slide_index(.x , .i = YearMonth, .f = ~expected_shortfall(., 
                                                                                                         "xret_index"), 
                                                                     .before = 2, # 24 months
                                                                     .after = -1/12, # 1 month
                                                                     .complete = TRUE))
  )
  ) %>% 
  select(-data) %>% 
  unnest(ES_list) %>% 
  group_by(FundId) %>% 
  nest() %>% 
  mutate(ES_list = map(data, ~ .x %>% 
                         mutate(ES_fund = slider::slide_index(.x , .i = YearMonth, .f = ~expected_shortfall(., 
                                                                                                              "xret"), 
                                                                .before = 2, # 24 months
                                                                .after = -1/12, # 1 month
                                                                .complete = TRUE))
  )
  ) %>% 
  select(-data) %>% 
  unnest(ES_list) %>%
  ungroup()


# Save file

save(fund_data_3, file = "pipeline/calculate_TailRisk/TailSens_and_ES.Rdata")  

# Now calculate TailRisk

load("pipeline/calculate_TailRisk/TailSens_and_ES.Rdata")

# Change NULL values to NA

fund_data_3[["TailSens"]][sapply(fund_data_3[["TailSens"]], is.null)] <- NA # Set NULL values in list to NA
fund_data_3[["ES_market"]][sapply(fund_data_3[["ES_market"]], is.null)] <- NA 
fund_data_3[["ES_fund"]][sapply(fund_data_3[["ES_fund"]], is.null)] <- NA 



fund_data_3 <- fund_data_3 %>% 
  mutate(TailSens = as.numeric(unlist(TailSens)),
         ES_market = as.numeric(unlist(ES_market)),
         ES_fund = as.numeric(unlist(ES_fund)),
         TailRisk = TailSens*(abs(ES_fund)/abs(ES_market)))

# Save file

df <- fund_data_3

save(df, file = "pipeline/calculate_TailRisk/TailRisk.Rdata")
