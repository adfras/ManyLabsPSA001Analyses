data {
  int<lower=1> N;            // observations
  int<lower=1> J;            // persons
  int<lower=1> K;            // sites
  int<lower=1> P;            // predictors (including intercept)
  matrix[N, P] X;
  array[N] int<lower=1,upper=J> person;
  array[N] int<lower=1,upper=K> site;
  // site membership per person (first site observed for that person)
  array[J] int<lower=1,upper=K> person_site;
  vector[N] y;
  // observation weights (1 = include in likelihood; 0 = hold out for CV)
  array[N] int<lower=0,upper=1> w;

  // holdout indices for this fold (so we only output log_lik for the holdout)
  int<lower=1> N_holdout;
  array[N_holdout] int<lower=1,upper=N> holdout_idx;
}

parameters {
  vector[P] beta;                    // global fixed effects

  vector<lower=0, upper=2>[P] tau_site;       // site offsets (shrunk)
  matrix[K, P] site_raw;

  vector<lower=0, upper=2>[P] tau_person;     // person offsets (shrunk)
  matrix[J, P] person_raw;

  real<lower=0.01, upper=5> sigma_base;       // base residual SD
  real<lower=0, upper=1> tau_sigma;           // person SD inflation
  vector[J] z_sigma;                 // person SD offsets
}

transformed parameters {
  matrix[K, P] site_eff;
  matrix[J, P] person_eff;
  vector<lower=0>[J] sigma_person;

  row_vector[P] tau_site_row = to_row_vector(tau_site);
  row_vector[P] tau_person_row = to_row_vector(tau_person);
  for (k in 1:K) site_eff[k] = site_raw[k] .* tau_site_row;
  for (j in 1:J) person_eff[j] = person_raw[j] .* tau_person_row;
  sigma_person = sigma_base * exp(z_sigma * tau_sigma);
}

model {
  // priors
  beta ~ normal(0, 1);
  tau_site ~ normal(0, 0.2);
  tau_person ~ normal(0, 0.3);
  to_vector(site_raw) ~ normal(0, 1);
  to_vector(person_raw) ~ normal(0, 1);

  sigma_base ~ lognormal(log(0.3), 0.25);
  tau_sigma ~ normal(0, 0.2);
  z_sigma ~ normal(0, 1);

  // likelihood (fit on training observations only)
  {
    vector[N] mu = X * beta;
    mu += rows_dot_product(X, site_eff[site]);
    mu += rows_dot_product(X, person_eff[person]);
    for (n in 1:N) {
      if (w[n] == 1) {
        target += normal_lpdf(y[n] | mu[n], sigma_person[person[n]]);
      }
    }
  }
}

generated quantities {
  // Holdout-only log-likelihood to keep CmdStan CSV outputs small.
  vector[N_holdout] log_lik_holdout;
  {
    vector[N] mu = X * beta;
    mu += rows_dot_product(X, site_eff[site]);
    mu += rows_dot_product(X, person_eff[person]);
    for (h in 1:N_holdout) {
      int n = holdout_idx[h];
      log_lik_holdout[h] = normal_lpdf(y[n] | mu[n], sigma_person[person[n]]);
    }
  }
}

