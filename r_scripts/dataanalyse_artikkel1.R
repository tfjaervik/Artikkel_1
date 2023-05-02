# Load packages -----------------------------------------------------------

library(zoo)
library(copula)
library(tidyverse)
library(purrr)
library(furrr)


# Load data ---------------------------------------------------------------

# Choose number of funds, and hence number of cores to use

numb_funds <- parallel::detectCores()*2 

# load("funds_cleaned.Rdata")

# load("data/stock_data_OBX.Rdata")
load("data/rdata/funds_list")

  # funds <- df3
# rm(df3)
# Helper function

`%nin%` <- function(x, table){ 
  match(x, table, nomatch = 0L) == 0L
}

# Extract names of funds that have been run

covered_funds <- 
  fs::dir_ls("data/ltd_estimates", type = "file") %>% 
  #fs::dir_ls("data/test_data", type = "file") %>% 
  fs::path_file() %>% 
  fs::path_ext_remove() %>% 
  str_remove("ltd_estimates_")


# Change name to funds_list and remove funds that have been run 

# funds_list <- mod_data_list
# rm(mod_data_list)



if(length(which(names(funds_list) %in% covered_funds)) > 0){
  funds_list <- funds_list[- which(names(funds_list) %in% covered_funds)]
}


# funds <- funds %>% 
#   filter(SecurityId %nin% covered_funds)





# Split the dataset on SecurityId and put into list

# funds_list <- split(funds, f = funds$SecurityId)


# Create example list

# funds_list_initial <- funds_list[1:64]
# 
# # Shorten length of dataframes in example list to 104
# 
# funds_list_initial <- funds_list_initial %>% 
#   map(~.x %>% slice_head(n = 106))

funds_list_initial <- funds_list

# Remove funds_list to free memory

rm(funds_list)


# File names  ---------------------------------------------

output_paths_initial <- 
  names(funds_list_initial) %>% 
  paste0("ltd_estimates_", .) %>% 
  fs::path("data", "ltd_estimates", 
           ., ext = "rds")
# fs::path("data", "ltd_estimates", ., ext = "rds")




# Lists needed  ---------------------------------------------

# Copula combinations

all_cop_types <- 
  list(LTD = c("clayton", "rot-gumbel", "rot-joe", "rot-galambos"),
       NTD = c("gauss", "frank", "fgm", "plackett"),
       UTD = c("gumbel", "joe", "galambos", "rot-clayton")) %>% 
  purrr::cross(.)

# Store optimized values

results_optim <- vector(mode = "list", length = length(all_cop_types))
names(results_optim) <- paste0("Combination", 1:length(all_cop_types))

# Store least nll (robustness test for copula selection)

results_optim_nll <- vector(mode = "list", length = length(all_cop_types))
names(results_optim_nll) <- paste0("Combination", 1:length(all_cop_types))

mix.cop <- vector(mode = "list", length = length(all_cop_types))
names(mix.cop) <- paste0("Opt", 1:length(all_cop_types))


# To calculate optimal mixture of copulas

copula <- list(weights = NULL, cop.objects = NULL, cop.names = NULL)

# List of functions to send to nlminb

# function_list <- sapply(1:length(all_cop_types), function(i)  
#   function(p, ...){
#     L(..., p, cop_type_LTD = all_cop_types[[i]]$LTD,
#       cop_type_NTD = all_cop_types[[i]]$NTD,
#       cop_type_UTD = all_cop_types[[i]]$UTD)
#   })

# Originally had "x" instead of "..." in L() above

# Note that the functions all seemingly depend on some unknown parameter i, 
# but the following code confirms that it works as intended

# as.list(environment(function_list[[2]]))
# as.list(environment(function_list[[3]]))



# Define functions --------------------------------------------------------
# Optimized copula densities ####

LTD_density_optimized <- function(x, t, cop_type_LTD = c("clayton", "rot-gumbel", 
                                                         "rot-joe", "rot-galambos")) {
  
  if(cop_type_LTD == "clayton"){
    
    return(rvinecopulib::dbicop(x, cop_type_LTD, parameters = t))
    
  } else if(cop_type_LTD == "rot-gumbel"){
    
    return(rvinecopulib::dbicop(x, sub(".*-", "", cop_type_LTD), parameters = t, rotation = 180))
    
  } else if(cop_type_LTD == "rot-joe"){
    
    return(rvinecopulib::dbicop(x, sub(".*-", "", cop_type_LTD), parameters = t, rotation = 180))
    
  } else if(cop_type_LTD == "rot-galambos"){
    
    return(dCopula(1-x, galambosCopula(t)))
    
  } else{
    
    stop("cop_type_LTD must be in c('clayton', 'rot-gumbel', 'rot-joe', 'rot-galambos')")
    
  }
}

NTD_density_optimized <- function(x, t, cop_type_NTD = c("gauss", "frank", 
                                                         "fgm", "plackett")){
  
  if(cop_type_NTD == "gauss"){
    
    return(rvinecopulib::dbicop(x, "gaussian", parameters = t))
    
  } else if(cop_type_NTD == "frank"){
    
    return(rvinecopulib::dbicop(x, "frank", parameters = t))
    
  } else if(cop_type_NTD == "fgm"){
    
    return(dCopula(x, fgmCopula(t)))
    
  } else if(cop_type_NTD == "plackett"){
    
    return(dCopula(x, plackettCopula(t)))
    
  } else{
    
    stop("cop_type_LTD must be in c('gauss', 'frank', 'fgm', 'plackett')")
    
  }
  
}


UTD_density_optimized <- function(x, t, cop_type_UTD = c("gumbel", "joe", 
                                                         "galambos", "rot-clayton")){
  
  if(cop_type_UTD == "gumbel"){
    
    return(rvinecopulib::dbicop(x, family = "gumbel", parameters = t))
    
  } else if(cop_type_UTD == "joe"){
    
    return(rvinecopulib::dbicop(x, family = "joe", parameters = t))
    
  } else if(cop_type_UTD == "galambos"){
    
    return(dCopula(x, galambosCopula(t)))
    
  } else if(cop_type_UTD == "rot-clayton"){
    
    return(rvinecopulib::dbicop(x, family = "clayton", rotation = 180,
                                parameters = t))
    
  } else{
    
    stop("cop_type_UTD must be in c('gumbel', 'joe', 'galambos', 'rot-clayton')")
    
  }
}

# Weights for optimization ####

# Used for intial guess of weights

delta.n2w <- function(m, delta){
  
  foo <- log(delta/delta[1])
  tdelta <- as.vector(tail(foo, m - 1))
  return(tdelta)
}

# Normalize weights so they sum to 1

delta.w2n <- function(m, tdelta){
  
  # set first element to one and fill in the last m - 1 elements with working parameters and take exp
  delta <- c(1, exp(tdelta))
  
  # normalize
  delta = delta/sum(delta)
  
  return(delta)
}


# Negative log-likelihood ####

# Optimized nll

L_optimized <- function(x, p, cop_type_LTD, cop_type_NTD, cop_type_UTD){
  
  
  # Transform copula parameters into valid values
  
  
  # LTD parameter transformation
  
  # Try a transformation from -1e-10 to 28 
  if(cop_type_LTD == "clayton"){
    coppar1 <- 1e-10 + 27.9*(exp(p[1])/(1 + exp(p[1])))
  } #Try a transformation from 1 to 50
  else if(cop_type_LTD == "rot-gumbel"){
    coppar1 <- 1 + 48.9*exp(p[1])/(1 + exp(p[1]))
  } # Try a transformation from 1 to 30
  else if(cop_type_LTD == "rot-joe"){
    coppar1 <- 1 + 28.9*exp(p[1])/(1 + exp(p[1]))
  }
  else if(cop_type_LTD == "rot-galambos"){
    coppar1 <- 50*(exp(p[1])/(1 + exp(p[1])))
  }
  
  # NTD parameter transformation
  
  if(cop_type_NTD == "gauss"){
    coppar2 <- sin(p[2])
  } # Try a transformation from -35 to 35
  else if(cop_type_NTD == "frank"){
    coppar2 <- -34.9 + 69.9*(exp(p[2])/(1 + exp(p[2])))
  }
  else if(cop_type_NTD == "fgm"){
    coppar2 <- sin(p[2]) # can only take values in [-1,1]
  }
  else if(cop_type_NTD == "plackett"){
    coppar2 <- 50*(exp(p[2])/(1 + exp(p[2]))) #from 0 to 50 to avoid NA evaluation
  }
  
  # UTD parameter transformation
  # Try a transformation from 1 to 50
  if(cop_type_UTD == "gumbel"){
    coppar3 <- 1 + 48.9*exp(p[3])/(1 + exp(p[3]))
  } # Try a transformation from 1 to 30
  else if(cop_type_UTD == "joe"){
    coppar3 <- 1 + 28.9*exp(p[3])/(1 + exp(p[3]))
  }
  else if(cop_type_UTD == "galambos"){
    coppar3 <- 50*(exp(p[3])/(1 + exp(p[3])))
  } # Try a transformation from -1e-10 to 28
  else if(cop_type_UTD == "rot-clayton"){
    coppar3 <- 1e-10 + 27.9*(exp(p[3])/(1 + exp(p[3])))
  }
  
  # Transform weights
  
  w <- delta.w2n(m = 3, tdelta = c(p[4], p[5]))
  
  w1 <- w[1]
  w2 <- w[2]
  w3 <- w[3]
  
  # Evaluate densities
  
  LTD <- LTD_density_optimized(x, t = coppar1, cop_type_LTD)
  
  NTD <- NTD_density_optimized(x, t = coppar2, cop_type_NTD)
  
  UTD <- UTD_density_optimized(x, t = coppar3, cop_type_UTD)
  
  # print(c(coppar1, coppar2, coppar3, w[1], w[2], w[3], sum(is.na(LTD)), sum(is.na(NTD)), sum(is.na(UTD))))
  # Define negative log-likelihood function
  
  nll <- -sum(log(w1*LTD + w2*NTD + w3*UTD + 1e-8)) # add 1e-8 to avoid log(0)
  
  return(nll)
}

# Initial parameter guess for optimization ####


init.guess <- function(x, cop_type_LTD, cop_type_NTD, cop_type_UTD){
  
  init <- rep(NA, 5)
  
  tau.n <- cor(x[,1], x[,2], method = "kendall")
  
  # Guess for first initial value
  
  if(cop_type_LTD == "clayton"){
    itau <- iTau(claytonCopula(), tau = tau.n)
    
    # Since the clayton copula parameter can take values in [-1, \infty), then
    # itau can take values in [-1, \infty). Since we are using log to (almost)
    # invert the parameter transformation in the nll function defined above, we 
    # must ensure a positive guess. Could also do it by writing
    # init[1] <- log(1 + itau)
    
    if(itau <= 0){
      init[1] <- log(0.1)
    }
    else {
      init[1] <- log(itau)
    }
  } 
  else if(cop_type_LTD == "rot-gumbel"){
    itau <- iTau(rotCopula(gumbelCopula()), tau = tau.n) # rotCopula seems to work here
    init[1] <- log(itau)
  }
  else if(cop_type_LTD == "rot-joe"){
    itau <- iTau(rotCopula(joeCopula()), tau = tau.n) # add more conditions as
    # in other script? # rotCopula seems to work here
    init[1] <- log(itau)
  }
  else if(cop_type_LTD == "rot-galambos"){
    itau <- iTau(rotCopula(galambosCopula()), tau = tau.n) # rotCopula seems to work here
    init[1] <- log(itau)
  }
  else {
    stop("not valid input for cop_type_LTD")
  }
  
  # Guess for second initial value
  
  if(cop_type_NTD == "gauss"){
    itau <- iTau(normalCopula(), tau = tau.n)
    # init[2] <- asin(itau) # Can go outside [-1,1], must fix
    init[2] <- itau # Try without transforming as -pi/2 <= asin(x) <= pi/2
    # the difference was negligible in simulation (equal to about fifth decimal)
  }
  else if(cop_type_NTD == "frank"){
    itau <- iTau(frankCopula(), tau = tau.n)
    if(itau <= 100){
      init[2] <- itau
    }
    else{
      init[2] <- 100
    }
  }
  
  
  # Note that for the fgm copula, Kendall's tau is restricted between -2/9 and
  # 2/9. Any values outside that range will be replaced automatically by the
  # corresponding lower/upper bound. See page 162 of An Introduction to Copulas
  # by Roger Nelsen (2006). R gives a warning for each time the fgm copula is 
  # included in the optimization (when tau.n is outside the range), 
  # but nothing goes wrong with the code.
  
  else if(cop_type_NTD == "fgm"){
    itau <- iTau(fgmCopula(), tau = tau.n)
    init[2] <- itau # Try without asin(x) transformation, as in the gaussian case
  }
  else if(cop_type_NTD == "plackett"){ # tau.n=1 | tau.n=-1 implies itau = NA
    itau <- iTau(plackettCopula(), tau = tau.n)
    if(all.equal(tau.n, 1) == TRUE){
      init[2] <- log(iTau(plackettCopula(), tau = 0.99))
    }
    else if(all.equal(tau.n, -1) == TRUE){
      log(iTau(plackettCopula(), tau = -0.99))
    }
    else{
      init[2] <- log(itau)
    }
    
  }
  else{
    stop("not valid input for cop_type_NTD")
  }
  
  # Guess for third initial value
  
  if(cop_type_UTD == "gumbel"){
    if(tau.n<= 0){
      itau <- iTau(gumbelCopula(), tau = 0)
      init[3] <- log(itau)
    }
    else{
      itau <- iTau(gumbelCopula(), tau = tau.n)
      init[3] <-  log(itau)
    }
  }
  
  # Notice that Joe density evaluates to inf/ NaN whenever copula parameter 
  # exceeds 4...(on dataset with two equal rows, i.e kendall's tau == 1)
  
  else if(cop_type_UTD == "joe"){ 
    itau <-iTau(joeCopula(), tau = tau.n)
    # if (itau <= 4){
    #   init[3] <- log(itau)
    # }
    # else{
    #   init[3] <- 4
    # }
    init[3] <- log(itau)
  }
  else if(cop_type_UTD == "galambos"){
    itau <- iTau(galambosCopula(), tau = tau.n)
    init[3] <- log(itau)
  }
  else if(cop_type_UTD == "rot-clayton"){
    itau <- iTau(rotCopula(claytonCopula()), tau = tau.n) # rotCopula seems to work here
    init[3] <- log(1 + itau) #ensures positve value in log(). See clayton above
  }
  else{
    stop("not valid input for cop_type_UTD")
  }
  
  # Initial guess for weights
  
  w <- delta.n2w(m = 3, rep(1/3, 3))
  #delta.w2n(m=3, c(0,0))
  
  init[4] <- w[1]
  init[5] <- w[2]
  
  return(init)
  
}

# Get results from copula parameter estimation using nlminb #####

# Get results optimized

get_results_optimized <- function(p, cop_type_LTD, cop_type_NTD, cop_type_UTD){
  
  vec <- rep(NA, 6)
  
  # NTD
  
  if(cop_type_LTD == "clayton"){
    theta_LTD <- 1e-10 + 27.9*(exp(p[1])/(1 + exp(p[1])))
  }
  else if(cop_type_LTD == "rot-gumbel"){
    theta_LTD <- 1 + 48.9*exp(p[1])/(1 + exp(p[1]))
  }
  else if(cop_type_LTD == "rot-joe"){
    theta_LTD <- 1 + 28.9*exp(p[3])/(1 + exp(p[3]))
  }
  else if(cop_type_LTD == "rot-galambos"){
    theta_LTD <- 50*(exp(p[1])/(1 + exp(p[1])))
  }
  else{
    stop("not valid input for cop_type_LTD")
  }
  
  vec[1] <- theta_LTD
  
  # NTD
  
  if(cop_type_NTD == "gauss"){
    theta_NTD <- sin(p[2])
  }
  else if(cop_type_NTD == "frank"){
    theta_NTD <- -34.9 + 69.9*(exp(p[2])/(1 + exp(p[2])))
  }
  else if(cop_type_NTD == "fgm"){
    theta_NTD <- sin(p[2])
  }
  else if(cop_type_NTD == "plackett"){
    theta_NTD <- 50*(exp(p[2])/(1 + exp(p[2])))
  }
  else{
    stop("not valid input for cop_type_NTD")
  }
  
  vec[2] <- theta_NTD
  
  # UTD
  
  if(cop_type_UTD == "gumbel"){
    theta_UTD <- 1 + 48.9*exp(p[3])/(1 + exp(p[3]))
  }
  else if(cop_type_UTD == "joe"){
    theta_UTD <- 1 + 28.9*exp(p[3])/(1 + exp(p[3]))
  }
  else if(cop_type_UTD == "galambos"){
    theta_UTD <- 50*(exp(p[3])/(1 + exp(p[3])))
  }
  else if(cop_type_UTD == "rot-clayton"){
    theta_UTD <- 1e-10 + 27.9*(exp(p[3])/(1 + exp(p[3])))
  }
  else{
    stop("not valid input for cop_type_UTD")
  }
  
  vec[3] <- theta_UTD
  
  # Weights
  
  w <- delta.w2n(m = 3, tdelta = c(p[4], p[5]))
  
  
  vec[4] <- w[1]
  vec[5] <- w[2]
  vec[6] <- w[3]
  
  # Result
  
  return(vec)
}

# Evaluate copula ####

# Evaluate copula optimized

evaluate_copula_optimized <- function(x, w,
                                      param_LTD, param_NTD, param_UTD,
                                      cop_name_LTD, cop_name_NTD, cop_name_UTD){
  
  # x is the data set
  # w is a weight vector
  # cop_name_xxx is just the names of the copulas
  # cop_type_xxx are copula objects. In this version of the function, they are
  # not only copula objects. For copulas in the package rvinecopulib this 
  # will be the parameter value of the corresponding copula instead of the copula
  # itself.
  
  # LTD evaluation
  
  if(cop_name_LTD == "clayton"){
    
    value1 <- rvinecopulib::pbicop(x, family = "clayton", parameters = param_LTD)
    
  } else if(cop_name_LTD == "rot-gumbel"){
    
    value1 <- rvinecopulib::pbicop(x, family = "gumbel", rotation = 180, 
                                   parameters = param_LTD)
    
  } else if(cop_name_LTD == "rot-joe"){
    
    value1 <- rvinecopulib::pbicop(x, family = "joe", rotation = 180,
                                   parameters = param_LTD)
  } else if(cop_name_LTD == "rot-galambos"){
    
    value1 <- (-1 + x[, 1] + x[, 2] + pCopula(1-x + 1e-8, 
                                              galambosCopula(param_LTD)))
    
  } else{
    
    stop("cop_name_LTD must be in c('clayton', 'rot-gumbel', 'rot-joe', 'rot-galambos')")
    
  }
  
  # NTD evaluation
  
  if(cop_name_NTD == "gauss"){
    
    value2 <- rvinecopulib::pbicop(x, family = "gaussian", 
                                   parameters = param_NTD)
    
  } else if(cop_name_NTD == "frank"){
    
    value2 <- rvinecopulib::pbicop(x, family = "frank", 
                                   parameters = param_NTD)
    
  } else if(cop_name_NTD == "fgm"){
    
    value2 <- pCopula(x, fgmCopula(param_NTD))
    
  } else if(cop_name_NTD == "plackett"){
    
    value2 <- pCopula(x, plackettCopula(param_NTD))
    
  } else{
    
    stop("cop_name_NTD must be in c('gauss', 'frank', 'fgm', 'plackett'")
    
  }
  
  
  # UTD evaluation
  
  if(cop_name_UTD == "gumbel"){
    
    value3 <- rvinecopulib::pbicop(x, family = "gumbel",
                                   parameters = param_UTD)
    
  } else if(cop_name_UTD == "joe"){
    
    value3 <- rvinecopulib::pbicop(x, family = "joe", 
                                   parameters = param_UTD)
    
  } else if(cop_name_UTD == "galambos"){
    
    value3 <- pCopula(x, param_UTD)
    
  } else if(cop_name_UTD == "rot-clayton"){
    
    value3 <- rvinecopulib::pbicop(x, family = "clayton", 
                                   rotation = 180, parameters = param_UTD)
    
  } else{
    
    stop("cop_name_UTD must be in c('gumbel', 'joe', 'galambos', rot-clayton')")
    
  }
  
  value <- w[1]*value1 + w[2]*value2 + w[3]*value3
  
  return(value)
  
}

# Calculate LTD ####


# Calculate LTD optimized

calc_LTD_optimized <- function(y){
  # Optimization ------------------------------------------------------------
  
  x <- y %>%
    #select(ret, ret_index) %>% # Select it in roll_calc_LTD below
    pobs() %>% # Remove this and it works fine.
    as.matrix()
  
  # x <- as.matrix(x)
  
  for (i in 1:length(all_cop_types)){
    
    l <- try(nlminb(init.guess(x, all_cop_types[[i]]$LTD, all_cop_types[[i]]$NTD,
                               all_cop_types[[i]]$UTD),
                    L_optimized,
                    x = x,
                    cop_type_LTD = all_cop_types[[i]]$LTD,
                    cop_type_NTD = all_cop_types[[i]]$NTD,
                    cop_type_UTD = all_cop_types[[i]]$UTD), silent = TRUE)
    
    
    
    if(class(l) == "try-error"){
      results_optim[[paste0("Combination", i)]] <- rep(NA_real_,6) # must be 
      # NA_real_ to work with the rest of the code, specifically calculation of 
      # the optimal copula mixtures. This is because copula objects must have 
      # parameters of class "numeric". Note that is.numeric(NA) = FALSE, while
      # is.na(NA_real_) = TRUE.
      
      results_optim_nll[[paste0("Combination", i)]] <- c(NA_real_)
      
    }else{
      results_optim[[paste0("Combination", i)]] <- get_results_optimized(l$par,
                                                                         
                                                                         all_cop_types[[i]]$LTD,
                                                                         all_cop_types[[i]]$NTD,
                                                                         all_cop_types[[i]]$UTD)
      
      
      results_optim_nll[[paste0("Combination", i)]] <- l$objective
      
    }
    
    
    names(results_optim[[i]]) <- c("coppar_LTD", "coppar_NTD", "coppar_UTD", 
                                   "weight_LTD", "weight_NTD", "weight_UTD")
    
    # Note that list element paste0("Combination", i) must match the naming given
    # in the "Lists needed" section of the script above. 
    
    
  }
  
  
  # Calculate optimal mixture of copulas ------------------------------------
  w.opt <- c() # Empty vector to fill with optimal weights in loop below
  
  for (i in 1:length(all_cop_types)){
    
    
    w.opt[1] <- results_optim[[i]]["weight_LTD"]
    w.opt[2] <- results_optim[[i]]["weight_NTD"]
    w.opt[3] <- results_optim[[i]]["weight_UTD"]
    
    names(w.opt) <- c("weight_LTD", "weight_NTD", "weight_UTD")
    
    # w.opt <- c(results_optim[[i]]["weight_LTD"] , 
    #            results_optim[[i]]["weight_NTD"],
    #            results_optim[[i]]["weight_UTD"])
    
    # LTD
    
    if (all_cop_types[[i]]$LTD == "clayton"){
      c1 <- claytonCopula(results_optim[[i]][["coppar_LTD"]])
      
    } else if (all_cop_types[[i]]$LTD == "rot-gumbel"){
      c1 <- gumbelCopula(results_optim[[i]][["coppar_LTD"]]) 
      
    } else if (all_cop_types[[i]]$LTD == "rot-joe"){
      c1 <- joeCopula(results_optim[[i]][["coppar_LTD"]]) 
      
    } else if (all_cop_types[[i]]$LTD == "rot-galambos"){
      c1 <- galambosCopula(results_optim[[i]][["coppar_LTD"]]) 
      
    }
    
    
    
    # NTD
    
    if (all_cop_types[[i]]$NTD == "gauss"){
      c2 <- normalCopula(results_optim[[i]][["coppar_NTD"]])
      
    } else if (all_cop_types[[i]]$NTD == "frank"){
      c2 <- frankCopula(results_optim[[i]][["coppar_NTD"]])
      
    } else if (all_cop_types[[i]]$NTD == "fgm"){
      c2 <- fgmCopula(results_optim[[i]][["coppar_NTD"]])
      
    } else if (all_cop_types[[i]]$NTD == "plackett"){
      c2 <- plackettCopula(results_optim[[i]][["coppar_NTD"]])
      
    }
    
    # UTD
    
    if (all_cop_types[[i]]$UTD == "gumbel"){
      c3 <- gumbelCopula(results_optim[[i]][["coppar_UTD"]])
      
    } else if (all_cop_types[[i]]$UTD == "joe"){
      c3 <- joeCopula(results_optim[[i]][["coppar_UTD"]])
      
    } else if (all_cop_types[[i]]$UTD == "galambos"){
      c3 <- galambosCopula(results_optim[[i]][["coppar_UTD"]])
      
    } else if(all_cop_types[[i]]$UTD == "rot-clayton"){
      c3 <- claytonCopula(results_optim[[i]][["coppar_UTD"]])
      
    }
    
    
    copula[["weights"]] <- w.opt
    copula[["cop.objects"]] <- c(LTD = c1, NTD = c2, UTD = c3)
    copula[["cop.names"]] <- c(LTD = all_cop_types[[i]]$LTD, 
                               NTD = all_cop_types[[i]]$NTD, 
                               UTD = all_cop_types[[i]]$UTD)
    
    
    mix.cop[[paste0("Opt", i)]] <- copula
    
  }
  
  
  
  
  
  # New way of calculating IAD ----------------------------------------------
  # Evaluation points
  
  vec1 <- seq(1:length(x[, 1]))/length(x[, 1])
  
  vec2 <- seq(1:length(x[, 2]))/length(x[, 2])
  
  ev_points <- expand.grid(vec1, vec2)
  ev_points <- as.matrix(ev_points)
  
  ec <- ec <- C.n(ev_points, X = x)
  ec <- ec[-length(ec)] # remove last point (must be same length as true below)
  # also, it is equal to the last element of true, which has to be removed, so 
  # removing them both do not add anything to the calculation.
  
  # Loop over estimated copulas and calculate IAD
  
  IAD <- rep(0, length = length(mix.cop))
  true <- rep(0, length = length(x[,1]))
  
  for (i in 1:length(mix.cop)){
    
    if(is.na(results_optim[[i]][["coppar_LTD"]]) == TRUE | 
       is.na(results_optim[[i]][["coppar_NTD"]]) == TRUE |
       is.na(results_optim[[i]][["coppar_UTD"]]) == TRUE |
       is.na(results_optim[[i]][["weight_LTD"]]) == TRUE |
       is.na(results_optim[[i]][["weight_NTD"]]) == TRUE |
       is.na(results_optim[[i]][["weight_UTD"]]) == TRUE )  {
      
      IAD[i] <- NA
      
    } else{
      
      
      # Evaluate copula j at the evaluation points
      
      # true <- try(evaluate.copula(ev_points, mix.cop[[i]][["weights"]], 
      #                             mix.cop[[i]][["cop.objects"]][["LTD"]], 
      #                             mix.cop[[i]][["cop.objects"]][["NTD"]],
      #                             mix.cop[[i]][["cop.objects"]][["UTD"]], 
      #                             mix.cop[[i]][["cop.names"]][["LTD"]], 
      #                             mix.cop[[i]][["cop.names"]][["UTD"]]),
      #             silent = TRUE)
      
      true <- try(evaluate_copula_optimized(ev_points, 
                                            c(results_optim[[i]]["weight_LTD"],
                                              results_optim[[i]]["weight_NTD"],
                                              results_optim[[i]]["weight_UTD"]), 
                                            results_optim[[i]][["coppar_LTD"]], 
                                            results_optim[[i]][["coppar_NTD"]],
                                            results_optim[[i]][["coppar_UTD"]], 
                                            all_cop_types[[i]]$LTD, 
                                            all_cop_types[[i]]$NTD,
                                            all_cop_types[[i]]$UTD),
                  silent = TRUE)
      
      if (class(true) == "try-error"){
        IAD[i] <- Inf
        
      } else {
        # Remove last element of true as it is equal to 1
        true <- true[-length(true)]
        
        IAD_aux <- ((ec-true)^2)/(true*(1-true))
        IAD[i] <- sum(IAD_aux)
      }
      
    }
  }
  
  # Calculate LTD -----------------------------------------------------------
  
  # Pick best model
  
  best.index <- which(rank(IAD) == 1)
  best.weights <- mix.cop[[best.index]][["weights"]]
  best.copula.LTD <- mix.cop[[best.index]][["cop.objects"]][["LTD"]]
  #best.copula.UTD <- mix.cop[[best.index]][["cop.objects"]][["UTD"]]
  
  
  # Calculate LTD of copulas in class LTD
  
  
  if (is.na(results_optim[[best.index]][["coppar_LTD"]]) == TRUE | 
      is.na(results_optim[[best.index]][["weight_LTD"]]) == TRUE){
    LTD <- NA
  } else if(startsWith(mix.cop[[best.index]][["cop.names"]][["LTD"]], "rot")){
    est.lambda.LTD <- lambda(mix.cop[[best.index]][["cop.objects"]][["LTD"]])[["upper"]]
    LTD <- best.weights[["weight_LTD"]]*est.lambda.LTD
  } else{
    est.lambda.LTD <- lambda(mix.cop[[best.index]][["cop.objects"]][["LTD"]])[["lower"]]
    LTD <- best.weights[["weight_LTD"]]*est.lambda.LTD
  }
  
  # Calculate UTD of copulas
  
  if (is.na(results_optim[[best.index]][["coppar_UTD"]]) == TRUE | 
      is.na(results_optim[[best.index]][["weight_UTD"]]) == TRUE){
    UTD <- NA
  } else if(startsWith(mix.cop[[best.index]][["cop.names"]][["UTD"]], "rot")){
    est.lambda.UTD <- lambda(mix.cop[[best.index]][["cop.objects"]][["UTD"]])[["lower"]]
    UTD <- best.weights[["weight_UTD"]]*est.lambda.UTD
  } else{
    est.lambda.UTD <- lambda(mix.cop[[best.index]][["cop.objects"]][["UTD"]])[["upper"]]
    UTD <- best.weights[["weight_UTD"]]*est.lambda.UTD
  }
  
  best.index.stored <- best.index
  # Must also calculate LTD based on copulas selected by the least nll
  # There is some copy and paste here from the code right above this line. 
  # Could make this code more robust by creating a function or mapping over a 
  # list
  
  # Best copula
  
  best.index.aux <- unlist(results_optim_nll)
  best.index <- which(rank(best.index.aux, ties.method = "random") == 1)
  
  # best.index <- which(unlist(results_optim_nll) == min(unlist(results_optim_nll)))
  best.weights <- mix.cop[[best.index]][["weights"]]
  best.copula.LTD <- mix.cop[[best.index]][["cop.objects"]][["LTD"]]
  #best.copula.UTD <- mix.cop[[best.index]][["cop.objects"]][["UTD"]]
  
  
  # Calculate LTD of copulas in class LTD
  
  
  if (is.na(results_optim[[best.index]][["coppar_LTD"]]) == TRUE |
      is.na(results_optim[[best.index]][["weight_LTD"]]) == TRUE){
    LTD_nll <- NA
  } else if(startsWith(mix.cop[[best.index]][["cop.names"]][["LTD"]], "rot")){
    est.lambda.LTD <- lambda(mix.cop[[best.index]][["cop.objects"]][["LTD"]])[["upper"]]
    LTD_nll <- best.weights[["weight_LTD"]]*est.lambda.LTD
  } else{
    est.lambda.LTD <- lambda(mix.cop[[best.index]][["cop.objects"]][["LTD"]])[["lower"]]
    LTD_nll <- best.weights[["weight_LTD"]]*est.lambda.LTD
  }
  
  return(paste0(LTD, "/", best.index.stored, "-", UTD, "-", LTD_nll, "-", best.index))
  
}



################################################################################
# Rolling calculation of LTD
################################################################################

# This function originally had rolling windows of 12 months. Changed it to 36
# months.

rolling_window_multivariate <- function(data, columns, function_name){
  
  # 'data' is the dataset. Must contain a column of type "Date"
  # 'columns' is a vector of column names, i.e. a character vector
  # 'function_name' is simply the name of the function (not with () after the name)
  # and not surrounded by "".
  
  
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
  ymu <- tail(unique(ym), -36)
  
  # Note: If there are not at least 12 months of data then ymu will be empty.
  # Hence, we remove any companies below that do not have at least 1 year of data.
  
  rowList <- lapply(ymu, function(x) which(ym > x-1 & ym < x))
  
  # Remove datasets with no rows in the row list
  
  keep_entries <- !(map_dbl(rowList, ~ length(.x)) == 0)
  
  rowList <- rowList[keep_entries]
  
  # Must remove the same elements from ymu
  
  ymu <- ymu[keep_entries]
  
  out <- data.frame(YearMonth = ymu,
                    n = lengths(rowList),
                    #monthly_estimate = sapply(rowList, aux_function))
                    monthly_estimate = map_chr(rowList, aux_function))
  
  data <- left_join(x = data, y = out, by = "YearMonth")
  
  return(data)
}

# Must rename column Date to date in each dataframe

df_list <- funds_list_initial %>% 
  map(., ~{ .x %>% 
      rename(date = Date)})


plan(multisession, workers = parallel::detectCores()*2)

# Start timing
ptm <- proc.time()


furrr::future_walk2(.x = df_list,
                    .y = output_paths_initial,
                    ~{
                      
                      rolling_window_multivariate(data = .x, 
                                                  columns = c("ret", "ret_index"),
                                                  function_name = calc_LTD_optimized) %>% 
                        readr::write_rds(file = .y)
                      
                    })

# Stop timing
proc.time() - ptm

plan(sequential)