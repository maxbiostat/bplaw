source("power_law_aux.r")
library(poweRlaw)
library(rstan)
rstan_options(auto_write = TRUE)
rstan:::rstudio_stanc("stan/continuous_power_law_diff_full.stan")
contplaw.fullDiff.model <- stan_model(file = "stan/continuous_power_law_diff_full.stan", allow_undefined = TRUE,
                               verbose = FALSE,
                               includes =  paste0('\n#include "',
                                                  file.path(getwd(), 'stan/diff_function_logintegral.hpp'), '"\n')
)