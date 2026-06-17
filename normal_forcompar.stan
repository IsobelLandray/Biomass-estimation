data {
  int<lower=1> K;               // number of bins (7)
  int<lower=0> y[K];            // observed counts in each bin
  int<lower=0> N;               // total count
  real bin_edges[K+1];          // bin boundaries
  real<lower=0> mean_prior;
}
parameters {
  real<lower=0> sigma;
  real <lower=18.5> mu;
  
}


transformed parameters {
  vector[K] p;
  vector[K] raw_p;
  
  for (k in 1:K) {
    real lower = bin_edges[k];
    real upper = bin_edges[k + 1];
    
      raw_p[k] = (normal_cdf(upper,mu,sigma) - normal_cdf(lower,mu,sigma));
  }
  
  p = raw_p / sum(raw_p);  
}




model {
  y ~ multinomial(p);           // likelihood
  sigma ~ normal(0, 7);         // prior for SD 
  mu ~ normal(mean_prior,5);     //prior for mean
}
