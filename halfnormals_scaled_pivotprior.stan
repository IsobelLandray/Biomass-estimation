data {
  int<lower=1> K;               // number of bins (7)
  int<lower=0> y[K];            // observed counts in each bin
  int<lower=0> N;               // total count
  real bin_edges[K+1];          // bin boundaries
  real<lower=0> pivot_prior1;
  real<lower=0> pivot_prior2;
}
parameters {
  real<lower=0> sigma1;
    real<lower=0> sigma2;
    real <lower=18.5> pivot;

}


transformed parameters {
  vector[K] p;
  vector[K] raw_p;

  for (k in 1:K) {
    real lower = bin_edges[k];
    real upper = bin_edges[k + 1];

    if (upper <= pivot) {
      // left side: include sigma1 factor
      raw_p[k] = 2 / (sigma1 + sigma2) * sigma1 * (normal_cdf((upper - pivot) / sigma1,0,1) - normal_cdf((lower - pivot) / sigma1,0,1));
    } else if (lower >= pivot) {
      // right side: include sigma2 factor
      raw_p[k] = 2 / (sigma1 + sigma2) * sigma2 *(normal_cdf((upper - pivot) / sigma2,0,1) - normal_cdf((lower - pivot) / sigma2,0,1));
    } else {
      // bin spans pivot: left portion then right portion (include their sigmas)
      raw_p[k] = 2 / (sigma1 + sigma2) * (sigma1 * (normal_cdf((pivot - pivot) / sigma1,0,1) - normal_cdf((lower - pivot) / sigma1,0,1)) + sigma2 * (normal_cdf((upper - pivot) / sigma2,0,1) - normal_cdf((pivot - pivot) / sigma2,0,1)));
    }
  }

  p = raw_p / sum(raw_p);  // optional if bin_edges cover all support; keeps numerical safety
}




model {
  y ~ multinomial(p);           // likelihood
  sigma1 ~ normal(0, 5);         // prior for SD for left end (i.e. smaller sd as steep drop)
  sigma2 ~ normal(0, 5);         // prior for SD for right end (i.e. larger sd as more wide)
  pivot ~ normal(pivot_prior1+pivot_prior2/2,5);
}

