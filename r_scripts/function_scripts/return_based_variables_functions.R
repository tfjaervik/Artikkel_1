################################################################################
#' Functions used to calculate return based variables
################################################################################



# Calculate CAPM-Alpha ----------------------------------------------------

# Note that mkt-rf from KF data library is calculated using value weighted 
# returns of all CRSP firms incorporated in the US and listed on the NYSE
# We have used S&P 500 returns as our proxy for market returns. To be consistent, 
# we thus continue using S&P 500 returns. Hence, we use xret_index  instead of
# mkt-rf.

# Recall that CAPM regression is xret = alpha + beta*xret_index

capm_alpha <- function(y){
  
  # y is the dataset in question
  # col1 is xret, col2 is xret_index
  
  y <- as_tibble(y)
  
  # aux <- lm(xret ~ xret_index, data = y)
  
  aux <- lm(xret ~ xret_index, data = y) # must later calculate with excess ret
  
  alpha <- coef(summary(aux))[1]
  
  return(alpha)
}

# Calculate FF3-Alpha ------------------------------------------------------

# Note that mkt-rf from KF data library is calculated using value weighted 
# returns of all CRSP firms incorporated in the US and listed on the NYSE
# We have used S&P 500 returns as our proxy for market returns. To be consistent, 
# we thus continue using S&P 500 returns. Hence, we use xret_index  instead of
# mkt-rf.

ff3_alpha <- function(y){
  
  y <- as_tibble(y)
  
  aux <- lm(xret ~ xret_index + smb + hml, data = y)
  
  alpha <- coef(summary(aux))[1]
  
  return(alpha)
}

# Calculate CAR-Alpha ---------------------------------------------------

car_alpha <- function(y){
  
  y <- as_tibble(y)
  
  aux <- lm(xret ~ xret_index + smb + hml + mom, data = y)
  
  alpha <- coef(summary(aux))[1]
  
  return(alpha)
}

# Calculate Beta ----------------------------------------------------------

calc_beta <- function(y){
  
  # first column is individual returns
  # second column is market return
  
  cov(y[, 1], y[, 2])/var(y[, 2])
}

# beta_df <- future_map(funds_list, 
#                       ~ {.x %>% 
#                           mutate(beta = rollapplyr(cbind(.x[["ret"]],
#                                                          .x[["ret_index"]]),
#                                                           104, calc_beta, fill = NA,
#                                                           by.column = FALSE, by = 4)) %>% 
#     fill(beta)
#   
# })




# Calculate downside beta -------------------------------------------------

calc_downside_beta <- function(y){
  
  # y is the dataset
  
  # This function is fed into a furrr loop which only operates on two columns, 
  # the first of which is individual returns and the second of which is market
  # returns. 
  
  mu <- mean(y[[2]]) # mean of market return
  
  conditional_data <- subset(y, y[[2]] < mu)
  
  # conditional_data <- y %>% 
  #   filter(y$ret_index < mu) # every dataset contains a column named ret-index
  
  return(calc_beta(conditional_data))
}



# downside_beta_df <- future_map(funds_list, 
#                       ~ {.x %>% 
#                           mutate(down_beta = rollapplyr(cbind(.x[["ret"]],
#                                                          .x[["ret_index"]]),
#                                                    104, calc_downside_beta, fill = NA,
#                                                    by.column = FALSE, by = 4)) %>% 
#                           fill(down_beta)
#                         
#                       })
# 





# Calculate coskewness ----------------------------------------------------

calc_coskewness <- function(y){
  
  # y is a matrix where the first column contains individual returns and the 
  # second column contain market returns
  
  mu_i <- mean(y[[1]])
  mu_m <- mean(y[[2]])
  
  coskew <- mean((y[[1]] - mu_i)*((y[[2]] - mu_m)^2)) / (sqrt(var(y[[1]]))*var(y[[2]]))
  
  return(coskew)
}

# coskewness_df <- future_map(funds_list_test, 
#                         ~{.x %>% 
#                            mutate(coskewness = rollapplyr(cbind(.x[["ret"]],
#                                                        .x[["ret_index"]]),
#                                                  104, calc_coskewness, fill = NA,
#                                                  by.column = FALSE, by = 4)) %>% 
#                     fill(coskewness)
# 
# })


# Calculate cokurtosis ----------------------------------------------------

calc_cokurtosis <- function(y){
  
  # y is a matrix where the first column contains individual returns and the 
  # second column contain market returns
  
  mu_i <- mean(y[[1]])
  mu_m <- mean(y[[2]])
  
  # mu_i <- mean(y[, 1])
  # mu_m <- mean(y[, 2])
  
  cokurt <- mean((y[[1]] - mu_i)*((y[[2]] - mu_m)^3)) / sqrt(var(y[[1]])*(var(y[[2]])^(3/2)))
  # cokurt <- mean((y[,1] - mu_i)*((y[, 2] - mu_m)^3)) / (sqrt(var(y[, 1]))*(var(y[, 2])^(3/2)))
  
  return(cokurt)
}

# cokurtosis_df <- future_map(funds_list_test, 
#                             ~{.x %>% 
#                                 mutate(cokurtosis = rollapplyr(cbind(.x[["ret"]],
#                                                                      .x[["ret_index"]]),
#                                                                104, calc_cokurtosis, fill = NA,
#                                                                by.column = FALSE, by = 4)) %>% 
#                                 fill(cokurtosis)
#                               
#                             })





# Calculate upside beta ---------------------------------------------------
calc_upside_beta <- function(y){
  
  # y is the dataset
  
  # This function is fed into a furrr loop which only operates on two columns, 
  # the first of which is individual returns and the second of which is market
  # returns. 
  
  mu <- mean(y[[2]]) # mean of market return
  
  conditional_data <- subset(y, y[[2]] > mu)
  
  # conditional_data <- y %>% 
  #   filter(y$ret_index < mu) # every dataset contains a column named ret-index
  
  return(calc_beta(conditional_data))
}


# 
# upside_beta_df <- future_map(funds_list_test, 
#                                ~ {.x %>% 
#                                    mutate(up_beta = rollapplyr(cbind(.x[["ret"]],
#                                                                        .x[["ret_index"]]),
#                                                                  104, calc_upside_beta, fill = NA,
#                                                                  by.column = FALSE, by = 4)) %>% 
#                                    fill(up_beta)
#                                  
#                                })



# data <- data_list[[1]] %>% head(n=15)
# columns <- c("xret", "xret_index", "hml", "smb")
# function_name <- ff3_alpha
# output_variable <- "ff3_alpha"


rolling_window_multivariate_modified <- function(data, columns, function_name, output_variable){
  
  # 'data' is the dataset. Must contain a column of type "Date"
  # 'columns' is a vector of column names, i.e. a character vector
  # 'function_name' is simply the name of the function (not with () after the name)
  # and not surrounded by "".
  
  if(!class(output_variable) == class("string")){
    stop("'output_variable' must be a string")
  }
  
  if(!("data.frame" %in% class(data))){
    stop("'data' must be a data frame")
  }
  
  if(!("Date" %in% sapply(data, class) & "yearmon" %in% sapply(data, class))){
    stop("Data must contain a column of type 'Date' and a column of type 'yearmon'")
  }
  
  if(!"YearMonth" %in% colnames(data)){
    stop("column of type 'yearmon' must be named 'YearMonth'")
  }
  
  if(class(columns) != "character" | !(all(columns %in% colnames(data)))){
    stop("Must provide valid columns names")
  }
  
  aux_function <- function(ix) {
    #print(ix)
    function_name(data[ix, columns])}
  
  #apply(data[ix, columns], 2, function_name)
  
  ym <- zoo::as.yearmon(data$date) 
  ymu <- tail(unique(ym), -12)
  
  rowList <- lapply(ymu, function(x) which(ym >= x-1 & ym < x))
  
  # If there are any missing years we have to remove them
  # Must not use which() here, because it may return an empty vector
  # a negative of which is also an empty vector. Instead, use logical operator
  # to denote which entries to keep
  
  keep_entries <- !(map_dbl(rowList, ~ length(.x)) == 0)
  
  rowList <- rowList[keep_entries]
  
  # Must remove same entry from ymu
  
  ymu <- ymu[keep_entries]
  
  out <- data.frame(YearMonth = ymu,
                    #n = lengths(rowList),
                    #monthly_estimate = sapply(rowList, aux_function))
                    monthly_estimate = map_chr(rowList, aux_function))
  
  
  names(out) <- c("YearMonth", output_variable)
  
  data <- left_join(x = data, y = out, by = "YearMonth")
  
  return(data)
}

