
functions{
    real power_law_lpdf(real x, real a, real m){
     return ( log(a-1)-log(m) -a *( log(x) - log(m)) );
    }
    /*TODO: implement rng for prior and posterior predictive checks*/
}
data{
  int<lower=0> K; // number of unique values  
  real values[K];
  int<lower=0> frequencies[K]; 
  real<lower=0> x_min;
  real<lower=0> alpha_shape;
  real<lower=0> alpha_rate;
}
parameters{
  real <lower=1> alpha;
}
model{
  /*WARNING: I know this distribution is CONTINUOUS and we shouldn't do this, CONCEPTUALLY*/
  /*HOWEVER: computationally, this makes a lot of sense. It is at least as efficient as not compressing*/
  for (k in 1:K) target += frequencies[k] * power_law_lpdf(values[k] | alpha, x_min);
  target += gamma_lpdf(alpha| alpha_shape, alpha_rate);
  // target += log(alpha - 1); // Jeffrey's prior
}