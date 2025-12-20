data {
  int<lower=1> N;
  int<lower=1> J;
  int<lower=1> K;
  int<lower=1> P;
  matrix[N, P] X;
  array[N] int<lower=1,upper=J> person;
  array[N] int<lower=1,upper=K> site;
  array[J] int<lower=1,upper=K> person_site;
  vector[N] y;
  // observation weights (1 = include in likelihood; 0 = hold out for CV)
  array[N] int<lower=0,upper=1> w;

  // holdout indices for this fold (so we only output log_lik for the holdout)
  int<lower=1> N_holdout;
  array[N_holdout] int<lower=1,upper=N> holdout_idx;
}

parameters {
  vector[P] beta;

  vector<lower=0, upper=2>[P] tau_site;
  matrix[K, P] site_raw;

  vector<lower=0, upper=2>[P] tau_person;
  matrix[J, P] person_raw;

  real<lower=0.01, upper=5> sigma;
}

transformed parameters {
  matrix[K, P] site_eff;
  matrix[J, P] person_eff;
  row_vector[P] tau_site_row = to_row_vector(tau_site);
  row_vector[P] tau_person_row = to_row_vector(tau_person);
  for (k in 1:K) site_eff[k] = site_raw[k] .* tau_site_row;
  for (j in 1:J) person_eff[j] = person_raw[j] .* tau_person_row;
}

model {
  beta ~ normal(0, 1);
  tau_site ~ normal(0, 0.2);
  tau_person ~ normal(0, 0.3);
  to_vector(site_raw) ~ normal(0, 1);
  to_vector(person_raw) ~ normal(0, 1);

  sigma ~ lognormal(log(0.3), 0.25);

  {
    vector[N] mu = X * beta;
    mu += rows_dot_product(X, site_eff[site]);
    mu += rows_dot_product(X, person_eff[person]);
    for (n in 1:N) {
      if (w[n] == 1) {
        target += normal_lpdf(y[n] | mu[n], sigma);
      }
    }
  }
}

generated quantities {
  vector[N_holdout] log_lik_holdout;
  {
    vector[N] mu = X * beta;
    mu += rows_dot_product(X, site_eff[site]);
    mu += rows_dot_product(X, person_eff[person]);
    for (h in 1:N_holdout) {
      int n = holdout_idx[h];
      log_lik_holdout[h] = normal_lpdf(y[n] | mu[n], sigma);
    }
  }
}

