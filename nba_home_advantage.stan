data {
  int<lower=1> G;                 // number of games
  int<lower=1> T;                 // number of teams
  int<lower=1> S;                 // number of seasons
  int<lower=0,upper=1> y[G];      // outcome: 1 = home win
  int<lower=1,upper=T> home_team[G];
  int<lower=1,upper=T> away_team[G];
  int<lower=1,upper=S> season[G];
  vector[G] zQ1;                  // standardized Q1 point diff
  vector[G] zQ2;                  // standardized Q2 point diff
}

parameters {
  real alpha0;                    // overall intercept (baseline home advantage)

  vector[S] alpha_raw;            // season deviations (standardized)
  real<lower=0> sigma_alpha;      // SD of season effects

  vector[T] theta;                // latent team strength
  real<lower=0> sigma_theta;      // SD of team strength

  vector[T] delta;                // team-specific home advantage
  real<lower=0> sigma_delta;      // SD of team home advantage

  real beta1;                     // effect of Q1 diff
  real beta2;                     // effect of Q2 diff
}

transformed parameters {
  vector[S] alpha;                // season intercepts
  alpha = alpha0 + sigma_alpha * alpha_raw;
}

model {
  // Hyperpriors
  alpha0      ~ normal(0, 2.5);
  sigma_alpha ~ normal(0, 1);     // half-normal via <0 constraint
  sigma_theta ~ normal(0, 1);
  sigma_delta ~ normal(0, 1);

  // Priors
  alpha_raw ~ normal(0, 1);
  theta     ~ normal(0, sigma_theta);
  delta     ~ normal(0, sigma_delta);

  beta1 ~ normal(0, 1);
  beta2 ~ normal(0, 1);

  // Likelihood
  for (g in 1:G) {
    real eta;
    eta = alpha[season[g]]
          + (theta[home_team[g]] - theta[away_team[g]])
          + delta[home_team[g]]
          + beta1 * zQ1[g]
          + beta2 * zQ2[g];

    y[g] ~ bernoulli_logit(eta);
  }
}

generated quantities {
  int y_rep[G];                   // posterior predictive outcomes

  for (g in 1:G) {
    real eta;
    eta = alpha[season[g]]
          + (theta[home_team[g]] - theta[away_team[g]])
          + delta[home_team[g]]
          + beta1 * zQ1[g]
          + beta2 * zQ2[g];

    y_rep[g] = bernoulli_logit_rng(eta);
  }
}
