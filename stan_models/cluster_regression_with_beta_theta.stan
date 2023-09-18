
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



