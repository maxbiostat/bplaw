source("power_law_aux.r")
library(poweRlaw)
library(rstan)
rstan_options(auto_write = TRUE)
# rstan:::rstudio_stanc("stan/continuous_power_law_diff.stan")
cplaw.model <- stan_model(file = "stan/continuous_power_law.stan")

