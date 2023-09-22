
data {
 int<lower = 0> N;
 int<lower = 0> N_cat1r;
 vector[N] y;
 int<lower = 1> cat1r[N];
}

parameters {

  real mu_base;
  real<lower=0> sigma_base; 
  real<lower=0> mu_diff;
  real<lower=0> sigma_diff;  
  real beta_base;

  real<lower=0> theta_cat1r_sd;
  real theta_cat1r_eps[N_cat1r];

}

transformed parameters {

  matrix[2, N] mu;
  vector[2] sigma;
  real<lower=0, upper=1> beta; 

  real cat1r_dev[N_cat1r]; 

  sigma[1] = sigma_base;
  sigma[2] = sigma_base + sigma_diff; 

  for (n in 1:N_cat1r) { 
   cat1r_dev[n] = theta_cat1r_sd * theta_cat1r_eps[n];
  }

  for (n in 1:N) {
   mu[1, n] = mu_base + cat1r_dev[cat1r[n]];
   mu[2, n] = mu[1, n] + mu_diff;
  }

  beta = inv_logit(beta_base); 

} 

model {

// --- Priors --- // 

 mu_base ~ normal(0, 2);
 mu_diff ~ normal(0, 2); 

 sigma_base ~ normal(0, 2);
 sigma_diff ~ normal(0, 2);

 beta_base ~ normal(0, 3);

 theta_cat1r_sd  ~ inv_gamma(8, 15);
 theta_cat1r_eps ~ normal(0, 3);


// --- Model --- // 

 for (n in 1:N)
   target += log_mix(beta,
                     normal_lpdf(y[n] | mu[2, n], sigma[2]),
                     normal_lpdf(y[n] | mu[1, n], sigma[1]));

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

   log_beta[1] += normal_lpdf(y[n] | mu[2, n], sigma[2]);
   log_beta[2] += normal_lpdf(y[n] | mu[1, n], sigma[1]);
   membership_l[, n] = exp(log_beta); 
 

  membership_p[1, n] = membership_l[1, n] / (membership_l[1, n] + membership_l[2, n]);
  membership_p[2, n] = membership_l[2, n] / (membership_l[1, n] + membership_l[2, n]);

  ind_sero[n] = binomial_rng(1, membership_p[1, n]);

  }

  pop_sero = sum(ind_sero);

}





