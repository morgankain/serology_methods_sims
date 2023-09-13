
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
  real theta_base;
  real theta_age;
  real beta_age;

}

transformed parameters {

  vector[2] mu;
  vector[2] sigma;
  vector[2] theta;

  mu[1] = mu_base;
  mu[2] = mu_base + mu_diff;
  sigma[1] = sigma_base;
  sigma[2] = sigma_base + sigma_diff; 

  theta[1] = inv_logit(theta_base);
  theta[2] = inv_logit(theta_base + theta_age);

} 

model {

 sigma_base ~ normal(0, 2);
 sigma_diff ~ normal(0, 2);
 mu_base ~ normal(0, 2);
 mu_diff ~ normal(0, 2); 
 theta_base ~ normal(0, 1);
 theta_age ~ normal(0, 1);

 beta_age ~ normal(0, 3);


 for (n in 1:N)
   target += log_mix(theta[age_index[n]],
                     normal_lpdf(y[n] | mu[1], sigma[1]),
                     normal_lpdf(y[n] | mu[2] + age[n] * beta_age, sigma[2]));


}


