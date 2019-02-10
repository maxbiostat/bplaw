library(rstan)
rstan_options(auto_write = TRUE)
rstan:::rstudio_stanc("stan/discrete_power_law.stan")
plaw.model <- stan_model(file = "stan/discrete_power_law.stan", allow_undefined = TRUE,
                         verbose =  TRUE,
                         includes = paste0('\n#include "',
                                           file.path(getwd(), 'stan/zeta.hpp'), '"\n')
)