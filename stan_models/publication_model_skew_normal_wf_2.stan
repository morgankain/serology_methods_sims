
functions {

real skew_mean(real loc, real scal, real shap) {
  real a;        
  real mm; 
  
  a = shap / (sqrt(1 + shap^2)); 
  mm = loc + scal * a * sqrt(2/pi()); 

  return mm;
 }

}

data {

  int<lower = 0> N;
  vector[N] y; 
  vector[N] cat1f;
  vector[N] cat2f;
  vector[N] con1f;

  real beta_base_prior_m;
  real beta_base_prior_v;

  real mu_diff_prior_m;
  real mu_diff_prior_v;
  real mu_pos_prior_m;
  real mu_pos_prior_v;
  real sigma_base_prior_m;  real sigma_diff_prior_m;  real sigma_base_prior_v;  real sigma_diff_prior_v;

  real skew_pos_prior_m;
  real skew_pos_prior_v; 

}

parameters {

  real<lower=0> mu_diff;
  real<lower=0> mu_pos;
  real<lower=0> sigma_base; 
  real<lower=0> sigma_diff;

  real<upper=0> skew_pos;

  real beta_base;

  real beta_cat1f_delta;
  real beta_cat2f_delta;
  real beta_con1f_delta; 

}

transformed parameters {

  vector[2] mu;
  vector[2] sigma;
  vector[N] beta_vec;

  sigma[1] = sigma_base;
  sigma[2] = sigma_base * sigma_diff; 

  for (n in 1:N) {  
   beta_vec[n] = inv_logit(beta_base + beta_cat1f_delta * cat1f[n] + beta_cat2f_delta * cat2f[n] + beta_con1f_delta * con1f[n]);
  }

  mu[1] = skew_mean(mu_pos, sigma[2], skew_pos) - mu_diff;
  mu[2] = mu_pos;

} 

model {

// --- Priors --- // 

 mu_diff ~ normal(mu_diff_prior_m, mu_diff_prior_v);
 mu_pos  ~ normal(mu_pos_prior_m, mu_pos_prior_v); 

 sigma_base ~ normal(sigma_base_prior_m, sigma_base_prior_v);
 sigma_diff ~ normal(sigma_diff_prior_m, sigma_diff_prior_v);

 skew_pos ~ normal(skew_pos_prior_m, skew_pos_prior_v);

 beta_base ~ normal(beta_base_prior_m, beta_base_prior_v);
 beta_cat1f_delta ~ normal(0, 3);
 beta_cat2f_delta ~ normal(0, 3);
 beta_con1f_delta ~ normal(0, 3);

// --- Model --- // 

 for (n in 1:N) {
   target += log_mix(beta_vec[n],
                     skew_normal_lpdf(y[n] | mu[2], sigma[2], skew_pos),
                     normal_lpdf(y[n] | mu[1], sigma[1]));
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

   log_beta[1, n] += skew_normal_lpdf(y[n] | mu[2], sigma[2], skew_pos);
   log_beta[2, n] += normal_lpdf(y[n] | mu[1], sigma[1]);

   membership_l[1, n] = exp(log_beta[1, n]); 
   membership_l[2, n] = exp(log_beta[2, n]); 

   membership_p[1, n] = membership_l[1, n] / (membership_l[1, n] + membership_l[2, n]);
   membership_p[2, n] = membership_l[2, n] / (membership_l[1, n] + membership_l[2, n]);

   ind_sero[n] = binomial_rng(1, membership_p[1, n]);

  }

  pop_sero = sum(ind_sero);

}


