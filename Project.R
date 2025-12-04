rm(list = ls())
setwd("C:/Users/LENOVO/Desktop/STAT 656/NBA/csv")

library(dplyr)
library(readr)
library(rstan)
library(bayesplot)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Read data

dat_raw <- read_csv("game_gameinfo_gamesummary_linescore_filtered.csv")

# Keep regular seasons 2014–2015 onward
dat <- dat_raw %>%
  filter(season_id >= 22014)

# Drop games where early quarter scores are missing
dat <- dat %>%
  filter(
    !is.na(pts_qtr1_home),
    !is.na(pts_qtr2_home),
    !is.na(pts_qtr1_away),
    !is.na(pts_qtr2_away)
  )

# Outcome: home win (1) / loss (0)

dat <- dat %>%
  mutate(
    home_win = ifelse(wl_home == "W", 1L, 0L)
  )

# Early momentum covariates

dat <- dat %>%
  mutate(
    dQ1 = pts_qtr1_home - pts_qtr1_away,
    dQ2 = pts_qtr2_home - pts_qtr2_away
  )

# Standardize (β’s on SD scale and α with nice interpretation)
dat <- dat %>%
  mutate(
    zQ1 = as.numeric(scale(dQ1)),
    zQ2 = as.numeric(scale(dQ2))
  )

# Indices for teams and seasons

# Use abbreviations to identify teams
teams <- sort(unique(c(dat$team_abbreviation_home,
                       dat$team_abbreviation_away)))
T <- length(teams)

dat <- dat %>%
  mutate(
    home_team = match(team_abbreviation_home, teams),
    away_team = match(team_abbreviation_away, teams)
  )

# Season indices
seasons <- sort(unique(dat$season_id))
S <- length(seasons)

dat <- dat %>%
  mutate(
    season_idx = match(season_id, seasons)
  )

# Build Stan data list

stan_data <- list(
  G    = nrow(dat),                  # number of games
  T    = T,                          # number of teams
  S    = S,                          # number of seasons
  y    = dat$home_win,               # 0/1 outcome
  home_team = dat$home_team,
  away_team = dat$away_team,
  season    = dat$season_idx,
  zQ1       = dat$zQ1,
  zQ2       = dat$zQ2
)

str(stan_data)

fit <- stan(
  file   = "nba_home_advantage.stan",
  data   = stan_data,
  chains = 4,
  iter   = 2000,
  warmup = 1000,
  seed   = 1234,
  control = list(adapt_delta = 0.9, max_treedepth = 12)
)

print(fit, pars = c("alpha0", "sigma_alpha", "sigma_theta",
                    "sigma_delta", "beta1", "beta2"),
      probs = c(0.1, 0.5, 0.9))

# Check convergence for the random effects
print(fit, pars = c("alpha", "theta", "delta"), probs = c(0.1, 0.5, 0.9))

mcmc_trace(as.array(fit), pars = c("alpha0", "beta1", "beta2"))
mcmc_pairs(as.array(fit), pars = c("alpha0", "beta1", "beta2"))

post <- rstan::extract(fit)

# alpha has dimension [iter, S]
alpha_draws <- post$alpha
season_home_prob <- apply(plogis(alpha_draws), 2, mean)
season_home_ci   <- apply(plogis(alpha_draws), 2,
                          quantile, probs = c(0.1, 0.5, 0.9))

season_summary <- data.frame(
  season_id = seasons,
  mean_home_p = season_home_prob,
  p10 = season_home_ci[1,],
  p50 = season_home_ci[2,],
  p90 = season_home_ci[3,]
)

season_summary

theta_draws <- post$theta   # [iter, T]
delta_draws <- post$delta   # [iter, T]

team_strength <- data.frame(
  team = teams,
  mean_theta = apply(theta_draws, 2, mean),
  p10_theta   = apply(theta_draws, 2, quantile, probs = 0.1),
  p50_theta   = apply(theta_draws, 2, quantile, probs = 0.5),
  p90_theta   = apply(theta_draws, 2, quantile, probs = 0.9),
  
  mean_delta = apply(delta_draws, 2, mean),
  p10_delta   = apply(delta_draws, 2, quantile, probs = 0.1),
  p50_delta   = apply(delta_draws, 2, quantile, probs = 0.5),
  p90_delta   = apply(delta_draws, 2, quantile, probs = 0.9)
)

team_strength %>% arrange(desc(mean_theta))  # rank teams by overall strength

beta1_draws <- post$beta1
beta2_draws <- post$beta2

quantile(beta1_draws, c(0.1, 0.5, 0.9))
quantile(beta2_draws, c(0.1, 0.5, 0.9))

y_rep <- post$y_rep   # [iter, G]

# Observed proportion of home wins
mean_obs <- mean(stan_data$y)

# Distribution of replicated proportions
prop_rep <- rowMeans(y_rep)

hist(prop_rep, breaks = 30, main = "PPC: proportion of home wins",
     xlab = "Proportion of home wins in replicated data")
abline(v = mean_obs, col = "red", lwd = 2)
