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
}

parameters {
  vector[P] beta;

  vector<lower=0, upper=2>[P] tau_site;
  matrix[K, P] site_raw;

  real<lower=0, upper=2> tau_person_intercept;
  vector[J] person_raw_intercept;

  real<lower=0.01, upper=5> sigma;
}

transformed parameters {
  matrix[K, P] site_eff;
  vector[J] person_eff_intercept;
  row_vector[P] tau_site_row = to_row_vector(tau_site);
  for (k in 1:K) site_eff[k] = site_raw[k] .* tau_site_row;
  person_eff_intercept = person_raw_intercept * tau_person_intercept;
}

model {
  beta ~ normal(0, 1);
  tau_site ~ normal(0, 0.2);
  tau_person_intercept ~ normal(0, 0.3);
  to_vector(site_raw) ~ normal(0, 1);
  person_raw_intercept ~ normal(0, 1);

  sigma ~ lognormal(log(0.3), 0.25);

  {
    vector[N] mu = X * beta;
    mu += rows_dot_product(X, site_eff[site]);
    mu += person_eff_intercept[person];
    for (n in 1:N) {
      if (w[n] == 1) {
        target += normal_lpdf(y[n] | mu[n], sigma);
      }
    }
  }
}

generated quantities {
  vector[N] log_lik;
  matrix[J, P] beta_person;

  {
    vector[N] mu = X * beta;
    mu += rows_dot_product(X, site_eff[site]);
    mu += person_eff_intercept[person];
    for (n in 1:N) {
      log_lik[n] = normal_lpdf(y[n] | mu[n], sigma);
    }
  }

  for (j in 1:J) {
    beta_person[j, 1] = beta[1] + site_eff[person_site[j], 1] + person_eff_intercept[j];
    beta_person[j, 2] = beta[2] + site_eff[person_site[j], 2];
  }
}
