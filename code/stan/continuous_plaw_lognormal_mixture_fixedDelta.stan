functions{
  real truncated_lognormal_lpdf(real x, real mu, real sigma, real x_min){
       real lconst = log(2) - log(erfc( (log(x_min)-mu)/(sqrt(2)*sigma) )) ; 
       real ldens = lognormal_lpdf(x | mu, sigma) ; 
       return(ldens + lconst);
  }
  real power_law_lpdf(real x, real a, real m){
     return ( log(a-1)-log(m) -a *( log(x) - log(m)) );
    }
  real pl_ln_mix_lpdf(real x, real alpha, real x_min, real mu, real sigma, real delta){
    real l1 = truncated_lognormal_lpdf(x | mu, sigma, x_min);
    real l2 = power_law_lpdf(x | alpha, x_min);
    real ans = log_mix(delta, l1, l2);
    return(ans);
  }
}
data{
  int<lower=0> K; // number of unique values  
  real values[K];
  int<lower=0> frequencies[K]; 
  real<lower=0> x_min;
  real<lower=0> alpha_shape;
  real<lower=0> alpha_rate;
  real<lower=0> adelta;
  real<lower=0> bdelta;
}
parameters{
  real<lower=1, upper=8> alpha;
  real<lower=0,upper=1> delta;
  real mu;
  real<lower=0> sigma;
}
model{
  /*WARNING: I know this distribution is CONTINUOUS and we shouldn't do this, CONCEPTUALLY*/
  /*HOWEVER: computationally, this makes a lot of sense. It is at least as efficient as not compressing*/
  for (k in 1:K) target += frequencies[k] * pl_ln_mix_lpdf(values[k] | alpha, x_min, mu, sigma, delta);
  target += gamma_lpdf(alpha| alpha_shape, alpha_rate);
  target += normal_lpdf(sigma | 0.5, 1.0);
  target += normal_lpdf(mu | 0.0, 2.0);
  target += beta_lpdf(delta | adelta, bdelta);
  // target += log(alpha - 1); // Jeffrey's prior
}
