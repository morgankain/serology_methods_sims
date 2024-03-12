
data {

  int<lower = 0> N;
  vector[N] y; 
  vector[N] cat1f;
  vector[N] cat2f;
  vector[N] con1f;

  real mu_base_prior;
  real mu_diff_prior;  real sigma_base_prior;  real sigma_diff_prior;

  real max_mfi;

}

parameters {

  real mu_base;
  real<lower=0> sigma_base; 
  real<lower=0> mu_diff;
  real<lower=0> sigma_diff;  

  real beta_base;

  real beta_cat1f_delta;
  real beta_cat2f_delta;
  real beta_con1f_delta; 

}

transformed parameters {

  matrix[2, N] mu;
  vector[2] sigma;
  vector[N] beta_vec;

  sigma[1] = sigma_base;
  sigma[2] = sigma_base + sigma_diff; 

  for (n in 1:N) {  
   beta_vec[n] = inv_logit(beta_base + beta_cat1f_delta * cat1f[n] + beta_cat2f_delta * cat2f[n] + beta_con1f_delta * con1f[n]);
  }

  for (n in 1:N) {
   mu[1, n] = mu_base;
   mu[2, n] = mu[1, n] + mu_diff;
  }

} 

model {

// --- Priors --- // 

 mu_base ~ normal(0, mu_base_prior);
 mu_diff ~ normal(0, mu_diff_prior); 

 sigma_base ~ normal(0, sigma_base_prior);
 sigma_diff ~ normal(0, sigma_diff_prior);

 beta_base ~ normal(0, 3);
 beta_cat1f_delta ~ normal(0, 3);
 beta_cat2f_delta ~ normal(0, 3);
 beta_con1f_delta ~ normal(0, 3);

// --- Model --- // 

 for (n in 1:N) {
   target += log_mix(beta_vec[n],
                     lognormal_lpdf(y[n] | mu[2, n], sigma[2]) - lognormal_lccdf(max_mfi | mu[2, n], sigma[2]),
                     lognormal_lpdf(y[n] | mu[1, n], sigma[1]) - lognormal_lccdf(max_mfi | mu[1, n], sigma[1]));
 }

}


generated quantities {

  matrix[2, N] membership_l;
  matrix[2, N] membership_p;
  array[N] int ind_sero;
  int pop_sero;
  matrix[2, N] log_beta;

  for (n in 1:N) {
   vector[2] beta_n;
 
   log_beta[1, n] = log(beta_vec[n]); 
   log_beta[2, n] = log(1 - beta_vec[n]); 

   log_beta[1, n] += lognormal_lpdf(y[n] | mu[2, n], sigma[2]) - lognormal_lcdf(max_mfi | mu[2, n], sigma[2]);
   log_beta[2, n] += lognormal_lpdf(y[n] | mu[1, n], sigma[1]) - lognormal_lcdf(max_mfi | mu[1, n], sigma[1]);

   membership_l[1, n] = exp(log_beta[1, n]); 
   membership_l[2, n] = exp(log_beta[2, n]); 

   membership_p[1, n] = membership_l[1, n] / (membership_l[1, n] + membership_l[2, n]);
   membership_p[2, n] = membership_l[2, n] / (membership_l[1, n] + membership_l[2, n]);

   ind_sero[n] = binomial_rng(1, membership_p[1, n]);

  }

  pop_sero = sum(ind_sero);

}


