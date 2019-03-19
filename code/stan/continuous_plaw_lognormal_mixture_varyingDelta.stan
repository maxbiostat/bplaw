functions{
  real truncated_lognormal_lpdf(real x, real mu, real sigma, real x_min){
       real lconst = log(2) - log(erfc( (log(x_min)-mu)/(sqrt(2)*sigma) )) ; 
       real ldens = lognormal_lpdf(x | mu, sigma) ; 
       return(ldens + lconst);
  }
  real power_law_lpdf(real x, real a, real m){
     return ( log(a-1)-log(m) -a *( log(x) - log(m)) );
    }
  real pl_ln_mix_lpdf(real x, real alpha, real x_min, real mu, real sigma, real phi){
    real l1 = truncated_lognormal_lpdf(x | mu, sigma, x_min);
    real l2 = power_law_lpdf(x | alpha, x_min);
    real delta = exp(-phi * x);
    real ans = log_mix(delta, l1, l2);
    return(ans);
  }
  real density2integrate(real x, real xc, real[] theta, real[] x_r, int[] x_i){
    real alpha = theta[1];
    real x_min = theta[2];
    real mu = theta[3];
    real sigma = theta[4];
    real phi = theta[5];
    real ldens = pl_ln_mix_lpdf(x | alpha, x_min, mu, sigma, phi);
    return(exp(ldens));
  }
}
data{
  int<lower=0> K; // number of unique values  
  real values[K];
  int<lower=0> frequencies[K]; 
  real<lower=0> x_min;
  real<lower=0> alpha_shape;
  real<lower=0> alpha_rate;
  real<lower=0> phi_shape;
  real<lower=0> phi_rate;
}
transformed data {
  real x_r[2];
  int x_i[10];
}
parameters{
  real<lower=1, upper=8> alpha;
  real<lower=0> phi;
  real mu;
  real<lower=0> sigma;
}
transformed parameters{
  real lconst = integrate_1d(density2integrate, x_min, positive_infinity(),  {alpha, x_min, mu, sigma, phi}, x_r, x_i, 0.01);
}
model{
  /*WARNING: I know this distribution is CONTINUOUS and we shouldn't do this, CONCEPTUALLY*/
  /*HOWEVER: computationally, this makes a lot of sense. It is at least as efficient as not compressing*/
  for (k in 1:K) target += frequencies[k] * ( pl_ln_mix_lpdf(values[k] | alpha, x_min, mu, sigma, phi) - lconst) ;
  target += gamma_lpdf(alpha| alpha_shape, alpha_rate);
  target += normal_lpdf(sigma | 0.5, 1.0);
  target += normal_lpdf(mu | 0.0, 2.0);
  target += gamma_lpdf(phi | phi_shape, phi_rate);
  // target += log(alpha - 1); // Jeffrey's prior
}
