functions{
    real diff_function_logintegral(real alpha, real z, real w, real m);
    real power_law_diff_lpdf(real x, real a, real m, real[] phi){
      real ldens = -a * log(x) + log1m_exp(-phi[1] - phi[2] * (x-1) - phi[3] * pow(x-1, 2) );
      return (ldens);
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
  real<lower=1> alpha;
  real<lower=0> phi[3];
}
transformed parameters{
    real lconst = log_diff_exp( (1-alpha)*log(x_min)-log(alpha-1),
  -phi[1] + diff_function_logintegral(alpha,
                                        exp(-phi[2]), // z
                                        exp(-phi[3]), // w
                                        x_min) );
}
model{
  /* WARNING: I know this distribution is CONTINUOUS and we shouldn't do this, CONCEPTUALLY */
  /* HOWEVER: computationally, this makes a lot of sense. It is at least as efficient as not compressing */
   for (k in 1:K) target += frequencies[k] * (power_law_diff_lpdf(values[k] | alpha, x_min, phi)-lconst);
  /* Priors */
  target += gamma_lpdf(alpha| alpha_shape, alpha_rate); // original prior by C. Gillespie was U(1, 8)
  for (i in 1:3) target += exponential_lpdf(phi[i] | 1.0/3.0);
}
