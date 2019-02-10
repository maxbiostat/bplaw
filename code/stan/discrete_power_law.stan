/*Assumes x_min = 1; if x_min > 1, we'd need a correction for that truncation*/
functions{
    real zeta(real s);
    /*TODO: implement rng for prior and posterior predictive checks*/
}
data{
  int<lower=0> K; // number of unique values
  int values[K];
  int<lower=0> frequencies[K];
  real<lower=0> alpha_shape;
  real<lower=0> alpha_rate;
}
parameters{
  real <lower=1> alpha;
}
model{
  real constant = log(zeta(alpha));
  for (k in 1:K) {
    target += frequencies[k] * (-alpha * log(values[k]) - constant);
  }
  target += gamma_lpdf(alpha |alpha_shape, alpha_rate);
} 
