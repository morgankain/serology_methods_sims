
data {
 int<lower = 0> N;
 vector[N] y;
}

parameters {

  real mu_base;
  real<lower=0> sigma_base; 
  real<lower=0> mu_diff;
  real<lower=0> sigma_diff;  
  real<lower=0, upper=1> beta;
}

transformed parameters {

  vector[2] mu;
  vector[2] sigma;

  mu[1] = mu_base;
  mu[2] = mu_base + mu_diff;
  sigma[1] = sigma_base;
  sigma[2] = sigma_base + sigma_diff; 

} 

model {

// --- Priors --- // 

 mu_base ~ normal(0, 2);
 mu_diff ~ normal(0, 2); 

 sigma_base ~ normal(0, 2);
 sigma_diff ~ normal(0, 2);

 beta ~ beta(5, 5);


// --- Model --- // 

 for (n in 1:N)
   target += log_mix(beta,
                     normal_lpdf(y[n] | mu[2], sigma[2]),
                     normal_lpdf(y[n] | mu[1], sigma[1]));

}

generated quantities {

  matrix[2, N] membership_l;
  matrix[2, N] membership_p;
  int ind_sero[N];
  int pop_sero;

  for (n in 1:N) {

   vector[2] beta_n;
 
   beta_n[1] = beta;
   beta_n[2] = 1 - beta;
   vector[2] log_beta = log(beta_n); 

   log_beta[1] += normal_lpdf(y[n] | mu[2], sigma[2]);
   log_beta[2] += normal_lpdf(y[n] | mu[1], sigma[1]);
   membership_l[, n] = exp(log_beta); 
 
   membership_p[1, n] = membership_l[1, n] / (membership_l[1, n] + membership_l[2, n]);
   membership_p[2, n] = membership_l[2, n] / (membership_l[1, n] + membership_l[2, n]);

   ind_sero[n] = binomial_rng(1, membership_p[1, n]);

  }

   pop_sero = sum(ind_sero);
 
}

