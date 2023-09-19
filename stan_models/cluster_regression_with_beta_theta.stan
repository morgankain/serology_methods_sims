
data {

  int<lower = 0> N;
  vector[N] y; 
  vector[N] age;
  int age_index[N];

}

parameters {

  real mu_base;
  real<lower=0> sigma_base; 
  real<lower=0> mu_diff;
  real<lower=0> sigma_diff;  
  real beta_base;
  real beta_age_delta;
  real mu_theta_age; 

}

transformed parameters {

  vector[2] mu;
  vector[2] sigma;
  vector[2] beta_vec;

  mu[1] = mu_base;
  mu[2] = mu_base + mu_diff;
  sigma[1] = sigma_base;
  sigma[2] = sigma_base + sigma_diff; 

  beta_vec[1] = inv_logit(beta_base);
  beta_vec[2] = inv_logit(beta_base + beta_age_delta);

} 

model {

// --- Priors --- // 

 mu_base ~ normal(0, 2);
 mu_diff ~ normal(0, 2); 

 sigma_base ~ normal(0, 2);
 sigma_diff ~ normal(0, 2);

 mu_theta_age ~ normal(0, 1);

 beta_base ~ normal(0, 3);
 beta_age_delta ~ normal(0, 3);


// --- Model --- // 

 for (n in 1:N)
   target += log_mix(beta_vec[age_index[n]],
                     normal_lpdf(y[n] | mu[1], sigma[1]),
                     normal_lpdf(y[n] | mu[2] + age[n] * mu_theta_age, sigma[2]));


}


generated quantities {

  matrix[2, N] membership_l;
  matrix[2, N] membership_p;
  int ind_sero[N];
  int pop_sero;

  for (n in 1:N) {
   vector[2] beta_n;
 
   beta_n[1] = beta_vec[age_index[n]];
   beta_n[2] = 1 - beta_vec[age_index[n]];
   vector[2] log_beta = log(beta_n); 

   log_beta[1] += normal_lpdf(y[n] | mu[1], sigma[1]);
   log_beta[2] += normal_lpdf(y[n] | mu[2] + age[n] * mu_theta_age, sigma[2]);
   membership_l[, n] = exp(log_beta); 

   membership_p[1, n] = membership_l[1, n] / (membership_l[1, n] + membership_l[2, n]);
   membership_p[2, n] = membership_l[2, n] / (membership_l[1, n] + membership_l[2, n]);

   ind_sero[n] = binomial_rng(1, membership_p[2, n]);

  }

  pop_sero = sum(ind_sero);

}





