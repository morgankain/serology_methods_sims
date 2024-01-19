
data {

 int<lower = 0> N;
 int<lower = 0> N_plate;
 int<lower = 0> N_loc;
 int<lower = 0> N_cov;
 int<lower = 0> N_spec;

 array[N] real<lower = 0> y;
 array[N] int<lower = 1> plate;
 array[N] int<lower = 1> loc;
 array[N] int<lower = 1> spec;

 matrix[N, N_cov] mm;

}

parameters {

  real mu_base;
  real beta_base;
  real<lower=0> sigma_base; 

  vector[N_cov] beta_diff;
  real<lower=0> mu_diff;
  real<lower=0> sigma_diff;  

  real<lower=0> mu_plate_sd;
  vector[N_plate] mu_plate_eps;

  real<lower=0> mu_loc_sd;
  vector[N_loc] mu_loc_eps;

  real<lower=0> beta_loc_sd;
  vector[N_loc] beta_loc_eps;

}

transformed parameters {

  matrix[2, N] mu;
  vector[2] sigma;
  vector[N] beta_vec;

  vector[N_plate] mu_plate_dev; 
  vector[N_loc] mu_loc_dev;
  vector[N_loc] beta_loc_dev;

  sigma[1] = sigma_base;
  sigma[2] = sigma_base + sigma_diff; 

  for (n in 1:N_plate) { 
   mu_plate_dev[n] = mu_plate_sd * mu_plate_eps[n];
  }

  for (n in 1:N_loc) { 
   mu_loc_dev[n] = mu_loc_sd * mu_loc_eps[n];
   beta_loc_dev[n] = beta_loc_sd * beta_loc_eps[n];
  }

  for (n in 1:N) {
   mu[1, n] = mu_base + mu_plate_dev[plate[n]] + mu_loc_dev[loc[n]];
   mu[2, n] = mu[1, n] + mu_diff;
  }

  for (n in 1:N) {  
   beta_vec[n] = inv_logit(beta_base + mm[n, ] * beta_diff + beta_loc_dev[loc[n]]);
  }

} 

model {

// --- Priors --- // 

 mu_base ~ normal(0, 2);
 mu_diff ~ normal(0, 2); 

 sigma_base ~ normal(0, 2);
 sigma_diff ~ normal(0, 2);

 beta_base ~ normal(0, 3);
 beta_diff ~ normal(0, 2);

 mu_plate_sd ~ inv_gamma(8, 15);
 mu_loc_sd ~ inv_gamma(8, 15);
 beta_loc_sd ~ inv_gamma(8, 15);

 mu_plate_eps ~ normal(0, 3);
 mu_loc_eps ~ normal(0, 3);
 beta_loc_eps ~ normal(0, 3);


// --- Model --- // 

 for (n in 1:N)
   target += log_mix(beta_vec[n],
                     normal_lpdf(y[n] | mu[2, n], sigma[2]),
                     normal_lpdf(y[n] | mu[1, n], sigma[1]));

}

generated quantities {

 matrix[2, N] membership_l; matrix[2, N] membership_p; matrix[N_loc, N_spec] pop_sero = rep_matrix(0, N_loc, N_spec); array[N] int ind_sero; matrix[2, N] log_beta;
 vector[2] temp_beta_vec_1;
 vector[2] temp_beta_vec_2;
 real max_val = -500; 
 for (n in 1:N) {    log_beta[1, n] = log(beta_vec[n]);   log_beta[2, n] = log(1 - beta_vec[n]);     log_beta[1, n] += normal_lpdf(y[n] | mu[2, n], sigma[2]);  log_beta[2, n] += normal_lpdf(y[n] | mu[1, n], sigma[1]);

  temp_beta_vec_1[1] = log_beta[1, n];
  temp_beta_vec_1[2] = max_val;

  temp_beta_vec_2[1] = log_beta[2, n];
  temp_beta_vec_2[2] = max_val;
    membership_l[1, n] = exp(max(temp_beta_vec_1));   membership_l[2, n] = exp(max(temp_beta_vec_2));     membership_p[1, n] = membership_l[1, n] / (membership_l[1, n] + membership_l[2, n]);  membership_p[2, n] = membership_l[2, n] / (membership_l[1, n] + membership_l[2, n]);    ind_sero[n] = binomial_rng(1, membership_p[1, n]);

  pop_sero[loc[n], spec[n]] = pop_sero[loc[n], spec[n]] + ind_sero[n];   }

}

