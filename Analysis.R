library(dplyr)
library(lubridate)
library(rstan)
library(ggplot2)
library(bayesplot)
library(patchwork)
library(loo)


options(mc.cores = parallel::detectCores())


# replace with your working directory
setwd("C:/Michel/00_PhD/04_Fall_2025/Bayes/Project/Project/stat656-final-proj-main/stat656-final-proj-main")

# read in data
dataw <- read.csv('game_gameinfo_gamesummary_linescore_filtered.csv')

table(dataw[, "season_id"])

data_filtered = dataw


#############################################################
### Variables:
#############################################################


table(data_filtered[, "wl_home"])
table(data_filtered[, "season_id"])

# for theta and delta

table(data_filtered[, "team_id_home"])
table(data_filtered[, "team_id_away"])

## we have 4299 games, 9 seasons and 30 teams.

df <- data_filtered

# outcome: home win
df$y <- as.integer(df$wl_home == "W")

# season index
df$season_f <- factor(df$season_id)
season <- as.integer(df$season_f)
S <- nlevels(df$season_f)

# team index 
all_teams <- sort(unique(c(df$team_id_home, df$team_id_away)))
team_map  <- setNames(seq_along(all_teams), all_teams)

home <- unname(team_map[as.character(df$team_id_home)])
away <- unname(team_map[as.character(df$team_id_away)])
T <- length(all_teams)

# stan data
stan_data <- list(
  N = nrow(df),
  y = df$y,
  S = S,
  season = season,
  T = T,
  home = home,
  away = away
)

c(N = stan_data$N, S = stan_data$S, T = stan_data$T, home_win_rate = mean(df$y))

#############################################################
### Arma el DF con Variables - con covariables
#############################################################

# Si Q1=7, el local ganó el Q1 por 7. Si Q2=−5, perdió el Q2 por 5.

table(data_filtered[, 'season_type'])
table(data_filtered[, 'game_status_text'])


# 0) filter
df <- data_filtered %>%
  filter(season_type == "Regular Season",
         game_status_text == "Final") %>%
  mutate(
    y = as.integer(wl_home == "W"),
    Q1 = pts_qtr1_home - pts_qtr1_away,
    Q2 = pts_qtr2_home - pts_qtr2_away
  ) %>%
  # remove duplicates
  filter(!is.na(y), !is.na(Q1), !is.na(Q2)) %>%
  distinct(game_id, .keep_all = TRUE)


# indices consistentes para season y equipos
season_ids <- sort(unique(df$season_id))
team_ids   <- sort(unique(c(df$team_id_home, df$team_id_away)))

df <- df %>%
  mutate(
    season    = as.integer(factor(season_id, levels = season_ids)),
    home_team = as.integer(factor(team_id_home, levels = team_ids)),
    away_team = as.integer(factor(team_id_away, levels = team_ids))
  )

# Q1 and Q2 standardization
df <- df %>%
  mutate(
    zQ1 = as.numeric(scale(Q1)),
    zQ2 = as.numeric(scale(Q2))
  )


# stan data
stan_data <- list(
  G = nrow(df),
  T = length(team_ids),
  S = length(season_ids),
  y = df$y,
  home_team = df$home_team,
  away_team = df$away_team,
  season = df$season,
  zQ1 = df$zQ1,
  zQ2 = df$zQ2
)

table(df$season_id)
table(df$y)

#############################################################
### MOdel definition
#############################################################

logistic_0 <- "data {
  int<lower=1> G;
  int<lower=1> T;
  int<lower=1> S;

  array[G] int<lower=0,upper=1> y;
  array[G] int<lower=1,upper=T> home_team;
  array[G] int<lower=1,upper=T> away_team;
  array[G] int<lower=1,upper=S> season;

  vector[G] zQ1;
  vector[G] zQ2;
}

parameters {
  real alpha0; // baseline home advantage (overall)

  vector[S] alpha_raw; // season deviations (std)
  real<lower=0> sigma_alpha;

  vector[T] theta_raw; // team strength (std)
  real<lower=0> sigma_theta;

  vector[T] delta_raw; // team-specific home advantage (std)
  real<lower=0> sigma_delta;

  real beta1;
  real beta2;
}

transformed parameters {
  vector[S] alpha; // season intercepts
  vector[T] theta; // centered strengths
  vector[T] delta; // centered home-adv

  // Sum-to-zero / centering for identifiability
  alpha = alpha0 + sigma_alpha * (alpha_raw - mean(alpha_raw));
  theta =        sigma_theta * (theta_raw - mean(theta_raw));
  delta =        sigma_delta * (delta_raw - mean(delta_raw));
}

model {
  // Priors - regularzation
  alpha0      ~ normal(0, 2.5);

  sigma_alpha ~ normal(0, 1);
  sigma_theta ~ normal(0, 1);
  sigma_delta ~ normal(0, 1);

  alpha_raw ~ normal(0, 1);
  theta_raw ~ normal(0, 1);
  delta_raw ~ normal(0, 1);

  beta1 ~ normal(0, 1);
  beta2 ~ normal(0, 1);

  // Linear predictor
  vector[G] eta =
      alpha[season]
    + theta[home_team] - theta[away_team]
    + delta[home_team]
    + beta1 * zQ1
    + beta2 * zQ2;

  y ~ bernoulli_logit(eta);
}

generated quantities {
  array[G] int y_rep;
  vector[G] log_lik;

  for (g in 1:G) {
    real eta_g =
        alpha[season[g]]
      + theta[home_team[g]] - theta[away_team[g]]
      + delta[home_team[g]]
      + beta1 * zQ1[g]
      + beta2 * zQ2[g];

    y_rep[g]  = bernoulli_logit_rng(eta_g);
    log_lik[g]= bernoulli_logit_lpmf(y[g] | eta_g);
  }
}
"

# Model 1 in slides

sm_0 <- stan_model(model_code = logistic_0)
#saveRDS(sm_0, 'compiled_sm_0.rds')

fit_model0 <- sampling(
  object = sm_0,
  data = stan_data,
  seed = 9876,
  chains = 4,
  iter = 5000,     
  warmup = 1000,
  control = list(adapt_delta = 0.9, max_treedepth = 12)
)

############################################################
## Analysis 
############################################################

post <- rstan::extract(fit_model0)

inv_logit <- function(z) 1/(1+exp(-z))

# redefine based on post: I-terations, T-eams, S-easons
I <- length(post$alpha0)
T <- ncol(post$theta_raw)
S <- ncol(post$alpha_raw)

# theta (I * T) centrado + escalado 
theta <- (post$theta_raw - rowMeans(post$theta_raw)) *
  matrix(post$sigma_theta, nrow = I, ncol = T)

# delta (I * T): centrado + escalado 
delta <- (post$delta_raw - rowMeans(post$delta_raw)) *
  matrix(post$sigma_delta, nrow = I, ncol = T)

# alpha_s (I * S): alpha0 + sigma_alpha * (alpha_raw - mean(alpha_raw)) 
# This is the same as: alpha_s ~ N(alpha0, sigma_alpha)
alpha_s <- matrix(post$alpha0, nrow = I, ncol = S) +
  matrix(post$sigma_alpha, nrow = I, ncol = S) * (post$alpha_raw - rowMeans(post$alpha_raw))

# Probability of win home by season 
# this is: average teams: theta=0 and delta=0 (because they are centered)
p_by_season_draws <- inv_logit(alpha_s)

p_home_avg_by_season <- colMeans(p_by_season_draws)

dim(theta); dim(delta); dim(alpha_s); length(p_home_avg_by_season)

################# labels

season_levels <- sort(unique(df$season_id))  
all_teams <- sort(unique(c(df$team_id_home, df$team_id_away)))

summ_ci <- function(mat) {
  data.frame(
    mean = apply(mat, 2, mean),
    p50  = apply(mat, 2, median),
    p025 = apply(mat, 2, quantile, probs = 0.025),
    p975 = apply(mat, 2, quantile, probs = 0.975)
  )
}

############## 1) Summary by season: baseline of home win
season_tbl <- data.frame(
  season_id = season_levels,
  p_home_avg = p_home_avg_by_season
)

season_tbl

season_ci <- data.frame(
  season_id = season_levels,
  mean = colMeans(p_by_season_draws),
  p50  = apply(p_by_season_draws, 2, median),
  p025 = apply(p_by_season_draws, 2, quantile, probs = 0.025),
  p975 = apply(p_by_season_draws, 2, quantile, probs = 0.975)
)
season_ci

########## 2) Summary by team: 

theta_sum <- summ_ci(theta)
delta_sum <- summ_ci(delta)
home_strength_sum <- summ_ci(theta + delta)

team_tbl <- data.frame(
  team_id = all_teams,

  theta_p50 = theta_sum$p50,
  theta_lo  = theta_sum$p025,
  theta_hi  = theta_sum$p975,

  delta_p50 = delta_sum$p50,
  delta_lo  = delta_sum$p025,
  delta_hi  = delta_sum$p975,

  home_strength_p50 = home_strength_sum$p50,
  home_strength_lo  = home_strength_sum$p025,
  home_strength_hi  = home_strength_sum$p975
)

team_tbl


########## 3) map team_id and labels
abbr_map <- unique(rbind(
  df[, c("team_id_home", "team_abbreviation_home")],
  setNames(df[, c("team_id_away", "team_abbreviation_away")],
           c("team_id_home","team_abbreviation_home"))
))
names(abbr_map) <- c("team_id","abbr")

team_tbl <- merge(team_tbl, abbr_map, by = "team_id", all.x = TRUE)
team_tbl <- team_tbl[order(-team_tbl$home_strength_p50), ]

head(team_tbl[, c("abbr","team_id","theta_p50","delta_p50","home_strength_p50")])


## Delta summary

# team summary: median and 95% IC 
delta_sum <- data.frame(
  team_index = 1:T,
  delta_mean = apply(delta, 2, mean),
  delta_p50  = apply(delta, 2, median),
  delta_p025 = apply(delta, 2, quantile, probs = 0.025),
  delta_p975 = apply(delta, 2, quantile, probs = 0.975)
)

# 4) mapear índice -> team_id -> abreviación
#all_teams <- sort(unique(c(df$team_id_home, df$team_id_away)))

#abbr_map <- unique(rbind(
#  df[, c("team_id_home", "team_abbreviation_home")],
#  setNames(df[, c("team_id_away", "team_abbreviation_away")],
#           c("team_id_home","team_abbreviation_home"))
#))
#names(abbr_map) <- c("team_id","abbr")

delta_tbl <- delta_sum %>%
  mutate(team_id = all_teams[team_index]) %>%
  left_join(abbr_map, by = "team_id") %>%
  arrange(desc(delta_p50))

delta_tbl

############################################################
## Convergence 
############################################################

# Rhat 
rhat_all <- summary(fit_model0)$summary[, "Rhat"]
head(rhat_all)

# Rhat for slids
summary(fit_model0)$summary[
  c("alpha0","sigma_alpha","sigma_theta","sigma_delta","beta1","beta2"),
  "Rhat"
]

# Perfect!!

#     alpha0 sigma_alpha sigma_theta sigma_delta       beta1       beta2 
#  0.9999099   1.0003936   1.0002888   1.0008353   0.9998693   0.9998812 

sumkey <- summary(fit_model0)$summary[c("alpha0","sigma_alpha","sigma_theta","sigma_delta","beta1","beta2"),
                                      c("n_eff","Rhat")]
sumkey


#load("Trabajo251205")


############################################################
## Variances 
############################################################

# 1) hierarchichal SDs 
sd_tbl <- as.data.frame(summary(fit_model0)$summary[
  c("sigma_alpha","sigma_theta","sigma_delta"),
  c("mean","sd","2.5%","50%","97.5%","n_eff","Rhat")
])
sd_tbl$param <- rownames(sd_tbl)
sd_tbl <- sd_tbl %>%
  select(param, everything()) %>%
  mutate( var_mean = mean^2,
          var_50   = `50%`^2,
          var_2.5  = `2.5%`^2,
          var_97.5 = `97.5%`^2 )
sd_tbl




############################################################
## Graphs 
############################################################


plot_team_effect <- function(team_tbl, effect = c("delta","theta","home_strength"),
                             title_size = 16, text_size = 12) {
  effect <- match.arg(effect)

  cols <- switch(effect,
    delta = c(p50="delta_p50", lo="delta_lo", hi="delta_hi"),
    theta = c(p50="theta_p50", lo="theta_lo", hi="theta_hi"),
    home_strength = c(p50="home_strength_p50", lo="home_strength_lo", hi="home_strength_hi")
  )

  plot_df <- team_tbl %>%
    filter(!is.na(abbr)) %>%
    transmute(
      abbr,
      p50 = .data[[cols["p50"]]],
      lo  = .data[[cols["lo"]]],
      hi  = .data[[cols["hi"]]]
    ) %>%
    arrange(p50) %>%
    mutate(abbr = factor(abbr, levels = abbr))

  ggplot(plot_df, aes(x = abbr, y = p50)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2) +
    geom_point(size = 2) +
    labs(x = "Team", y = paste0(effect, " (log-odds)")) +
    theme_bw() +
    theme(
      axis.title.x = element_text(size = title_size),
      axis.title.y = element_text(size = title_size),
      axis.text.x  = element_text(size = text_size, angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y  = element_text(size = text_size)
    )
}

plot_team_effect(team_tbl, "delta", title_size = 18, text_size = 12)

plot_team_effect(team_tbl, "theta", title_size = 18, text_size = 13) + coord_flip()
plot_team_effect(team_tbl, "delta", title_size = 18, text_size = 13) + coord_flip()

plot_team_effect(team_tbl, "delta")
plot_team_effect(team_tbl, "theta")
plot_team_effect(team_tbl, "home_strength")


################## Betas


beta_df <- tibble(beta1 = post$beta1, beta2 = post$beta2) %>%
  pivot_longer(everything(), names_to = "param", values_to = "value") %>%
  mutate(param = factor(param, levels = c("beta1","beta2")))

ggplot(beta_df, aes(x = value, y = after_stat(density), color = param)) +
  geom_density(linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Coefficient (log-odds)", y = "Posterior density", color = NULL) +
  scale_color_discrete(
    labels = c(expression(beta[1]), expression(beta[2]))
  ) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text  = element_text(size = 12),
        legend.text = element_text(size = 14))

################## Alphas

alpha_s <- matrix(post$alpha0, nrow = I, ncol = S) +
  matrix(post$sigma_alpha, nrow = I, ncol = S) * (post$alpha_raw - rowMeans(post$alpha_raw))

season_levels <- sort(unique(df$season_id))

alpha_tbl <- data.frame(
  season = season_levels,
  p50  = apply(alpha_s, 2, median),
  lo   = apply(alpha_s, 2, quantile, probs = 0.025),
  hi   = apply(alpha_s, 2, quantile, probs = 0.975)
) %>%
  arrange(season) %>%
  mutate(season = factor(season, levels = season))

ggplot(alpha_tbl, aes(x = season, y = p50)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15) +
  geom_point(size = 2) +
  labs(x = "Season", y = expression(alpha[s]~"(log-odds)"),
       title = "") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text  = element_text(size = 12),
        plot.title = element_text(size = 16))

#League-wide home advantage by season (posterior median and 95% CrI)


####################################################################
####################################################################
####################################################################
## Predictive checks: this is goodnes of fit: is my model reproducing my data?
####################################################################
####################################################################
####################################################################

#post <- rstan::extract(fit_model0)
y    <- stan_data$y
yrep <- post$y_rep # iter * G

# 300 draws

set.seed(9876)
nd <- 300
idx <- sample(seq_len(nrow(yrep)), nd)
yrep_s <- yrep[idx, ]


bayesplot::ppc_stat(
  y = y,
  yrep = yrep_s,
  stat = "mean"
) + ggtitle("PPC: overall home-win rate")


bayesplot::ppc_stat(y, yrep_s, stat = "sd") + ggtitle("PPC: SD of outcomes")


season_g <- stan_data$season  

#S <- stan_data$S

# observed
obs_season <- sapply(1:S, function(s) mean(y[season_g == s]))

# replicated: nd x S
rep_season <- sapply(1:S, function(s) rowMeans(yrep_s[, season_g == s, drop = FALSE]))

season_ppc <- data.frame(
  season = 1:S,
  obs    = obs_season,
  lo     = apply(rep_season, 2, quantile, .025),
  mid    = apply(rep_season, 2, median),
  hi     = apply(rep_season, 2, quantile, .975)
)

ggplot(season_ppc, aes(x = season)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) +
  geom_line(aes(y = mid)) +
  geom_point(aes(y = obs), size = 2) +
  labs(x = "Season index", y = "Home-win rate", title = "PPC by season") +
  theme_bw()


home_g <- stan_data$home_team 



obs_home_team <- sapply(1:T, function(t) mean(y[home_g == t]))

rep_home_team <- sapply(1:T, function(t) rowMeans(yrep_s[, home_g == t, drop = FALSE]))

team_ppc <- data.frame(
  team = 1:T,
  obs  = obs_home_team,
  lo   = apply(rep_home_team, 2, quantile, .025),
  mid  = apply(rep_home_team, 2, median),
  hi   = apply(rep_home_team, 2, quantile, .975)
)

ggplot(team_ppc, aes(x = team)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) +
  geom_line(aes(y = mid)) +
  geom_point(aes(y = obs), size = 1.8) +
  labs(x = "Team index", y = "Home-win rate", title = "PPC by home team") +
  theme_bw()


# usar draws subsampleadas idx también
alpha <- post$alpha[idx, ]   # nd x S
theta <- post$theta[idx, ]   # nd x T
delta <- post$delta[idx, ]   # nd x T

beta1 <- post$beta1[idx]
beta2 <- post$beta2[idx]

zQ1 <- stan_data$zQ1
zQ2 <- stan_data$zQ2
away_g <- stan_data$away_team

eta_draws <- alpha[, season_g] +
  theta[, home_g] - theta[, away_g] +
  delta[, home_g] +
  beta1 %o% zQ1 +
  beta2 %o% zQ2

p_hat <- colMeans(plogis(eta_draws))  # posterior mean prob por juego

cal_df <- data.frame(y = y, p = p_hat) %>%
  mutate(bin = ntile(p, 10)) %>%
  group_by(bin) %>%
  summarise(p_mean = mean(p), y_mean = mean(y), n = n(), .groups = "drop")

ggplot(cal_df, aes(x = p_mean, y = y_mean)) +
  geom_point(size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "Mean predicted P(home win)", y = "Observed home-win rate",
       title = "Calibration (deciles)") +
  theme_bw()

brier_obs <- mean((y - p_hat)^2)

brier_rep <- apply(yrep_s, 1, function(yy) mean((yy - p_hat)^2))

hist(brier_rep, breaks = 30, main = "PPC: Brier score", xlab = "Brier(y_rep, p_hat)")
abline(v = brier_obs, lty = 2, lwd = 2)

###########################################################
### PPC agregado con log-likelihood total (o deviance)
###########################################################


# log-likelihood observado por draw

ll_draw <- rowSums(post$log_lik)

hist(ll_draw, breaks = 40,
     main = "Posterior distribution of total log-likelihood",
     xlab = "sum_g log p(y_g | params)")
abline(v = median(ll_draw), lty = 2)


## PPC agregado clásico 

y    <- stan_data$y
yrep <- post$y_rep

set.seed(9876)
idx <- sample(seq_len(nrow(yrep)), 500)
yrep_s <- yrep[idx, ]

T_obs <- mean((y - p_hat)^2)
T_rep <- apply(yrep_s, 1, function(yy) mean((yy - p_hat)^2))

hist(T_rep, breaks = 30, main = "PPC (global): Brier score",
     xlab = "T(y_rep)")
abline(v = T_obs, lty = 2, lwd = 2)


# Log score PPC 

eps <- 1e-9
T_obs <- -mean(y*log(p_hat+eps) + (1-y)*log(1-p_hat+eps))

T_rep <- apply(yrep_s, 1, function(yy) {
  -mean(yy*log(p_hat+eps) + (1-yy)*log(1-p_hat+eps))
})

hist(T_rep, breaks = 30, main = "PPC (global): log score",
     xlab = "T(y_rep)  (lower = better)")
abline(v = T_obs, lty = 2, lwd = 2)

# 3) un PPC global 

(ppp <- mean(T_rep >= T_obs))

####################################################################
## PPC: mean and variance (global, season, team)
####################################################################


set.seed(9876)
nd  <- 300
idx <- sample(seq_len(nrow(yrep)), nd)
yrep_s <- yrep[idx, ]

# -----------------------------
# 1) PPC GLOBAL: meand and variance
# -----------------------------
p_global_mean <- ppc_stat(y = y, yrep = yrep_s, stat = "mean") +
  ggtitle("PPC (global): mean(y) = overall home-win rate") +
  theme(legend.position = "none")

p_global_var <- ppc_stat(y = y, yrep = yrep_s, stat = var) +
  ggtitle("PPC (global): var(y)") +
  theme(legend.position = "none")

p_global_mean
p_global_var


# functions
mk_band_df <- function(rep_mat, obs_vec, group_ids, G_levels, label = "group") {
  obs_mean <- sapply(1:G_levels, function(k) mean(obs_vec[group_ids == k]))
  obs_var  <- sapply(1:G_levels, function(k) var (obs_vec[group_ids == k]))
  rep_mean <- sapply(1:G_levels, function(k) {
    rowMeans(rep_mat[, group_ids == k, drop = FALSE])
  })
  rep_var  <- sapply(1:G_levels, function(k) {
    apply(rep_mat[, group_ids == k, drop = FALSE], 1, var)
  })

  df_mean <- data.frame(
    group = 1:G_levels,
    obs   = obs_mean,
    lo    = apply(rep_mean, 2, quantile, .025),
    mid   = apply(rep_mean, 2, median),
    hi    = apply(rep_mean, 2, quantile, .975),
    stat  = "Mean"
  )

  df_var <- data.frame(
    group = 1:G_levels,
    obs   = obs_var,
    lo    = apply(rep_var, 2, quantile, .025),
    mid   = apply(rep_var, 2, median),
    hi    = apply(rep_var, 2, quantile, .975),
    stat  = "Variance"
  )

  out <- bind_rows(df_mean, df_var) %>%
    mutate(!!label := group) %>%
    select(-group)

  out
}

# -----------------------------
# 2) PPC BY SEASON: mean and variance
# -----------------------------

season_ppc <- mk_band_df(yrep_s, y, season_g, S, label = "season")

p_season_mean <- season_ppc %>% filter(stat == "Mean") %>%
  ggplot(aes(x = season)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) +
  geom_line(aes(y = mid)) +
  geom_point(aes(y = obs), size = 2) +
  labs(x = "Season", y = "Mean(y)", title = "PPC by season: mean(y)") +
  theme_bw()

p_season_var <- season_ppc %>% filter(stat == "Variance") %>%
  ggplot(aes(x = season)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) +
  geom_line(aes(y = mid)) +
  geom_point(aes(y = obs), size = 2) +
  labs(x = "Season", y = "Var(y)", title = "PPC by season: var(y)") +
  theme_bw()

p_season_mean
p_season_var


# -----------------------------
# 3) PPC BY TEAm (home_team): mean and variance
# -----------------------------

team_ppc <- mk_band_df(yrep_s, y, home_g, T, label = "team")

# Para slides: mejor ordenar equipos por el observado y usar pointrange + punto
p_team_mean <- team_ppc %>% filter(stat == "Mean") %>%
  mutate(team_f = reorder(factor(team), obs)) %>%
  ggplot(aes(x = team_f)) +
  geom_pointrange(aes(y = mid, ymin = lo, ymax = hi), linewidth = 0.3) +
  geom_point(aes(y = obs), size = 1.5, color = "red") +
  coord_flip() +
  labs(x = "Team (sorted) ", y = "Mean(y)",
       title = "PPC by home team (95% PI + median; dot = observed)") +
  theme_bw()

p_team_var <- team_ppc %>% filter(stat == "Variance") %>%
  mutate(team_f = reorder(factor(team), obs)) %>%
  ggplot(aes(x = team_f)) +
  geom_pointrange(aes(y = mid, ymin = lo, ymax = hi), linewidth = 0.3) +
  geom_point(aes(y = obs), size = 1.5, color = "red") +
  coord_flip() +
  labs(x = "Team (sorted)", y = "Var(y)",
       title = "PPC by home team (95% PI + median; dot = observed)") +
  theme_bw()

p_team_mean
p_team_var


### Junta todo

panel6 <-
  (p_global_mean | p_global_var) /
  (p_season_mean | p_season_var) /
  (p_team_mean   | p_team_var) +
  plot_layout(heights = c(1, 1, 2)) +
  plot_annotation(title = "")

panel6



##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
### Alternative model  - NO RUN, it is exploratory
##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################

logistic_1 <- "
data {
  int<lower=1> G;
  int<lower=1> T;
  int<lower=1> S;

  array[G] int<lower=0,upper=1> y;
  array[G] int<lower=1,upper=T> home_team;
  array[G] int<lower=1,upper=T> away_team;
  array[G] int<lower=1,upper=S> season;
}

parameters {
  real alpha0;

  vector[S] alpha_raw;
  real<lower=0> sigma_alpha;

  matrix[T,S] theta_raw;     // team strength by season (std)
  real<lower=0> sigma_theta; // global pooling
}

transformed parameters {
  vector[S] alpha;
  matrix[T,S] theta;         // centered within season

  alpha = alpha0 + sigma_alpha * (alpha_raw - mean(alpha_raw));

  for (s in 1:S) {
    vector[T] th_s = to_vector(theta_raw[, s]);
    theta[, s] = sigma_theta * (th_s - mean(th_s));
  }
}

model {
  // priors
  alpha0 ~ normal(0, 2.5);
  sigma_alpha ~ normal(0, 1);
  sigma_theta ~ normal(0, 1);

  alpha_raw ~ normal(0, 1);
  to_vector(theta_raw) ~ normal(0, 1);

  // likelihood
  for (g in 1:G) {
    real eta =
      alpha[season[g]]
      + theta[home_team[g], season[g]]
      - theta[away_team[g], season[g]];

    y[g] ~ bernoulli_logit(eta);
  }
}

generated quantities {
  array[G] int y_rep;
  vector[G] log_lik;

  for (g in 1:G) {
    real eta =
      alpha[season[g]]
      + theta[home_team[g], season[g]]
      - theta[away_team[g], season[g]];

    y_rep[g]  = bernoulli_logit_rng(eta);
    log_lik[g]= bernoulli_logit_lpmf(y[g] | eta);
  }
}
"


sm_1 <- stan_model(model_code = logistic_1)
 

fit_model1 <- sampling(
  object = sm_1,
  data = stan_data,
  seed = 9876,
  chains = 4,
  iter = 5000,     
  warmup = 1000,
  control = list(adapt_delta = 0.9, max_treedepth = 12)
)


print(fit_model1, pars=c("alpha0","sigma_alpha","sigma_theta"), probs=c(.025,.5,.975))
summary(fit_model1)$summary[c("alpha0","sigma_alpha","sigma_theta"), c("mean","sd","n_eff","Rhat")]


# Draws of theta (iter * T * S)
post1 <- rstan::extract(fit_model1, pars = "theta")
theta_draws <- post1$theta
I <- dim(theta_draws)[1]
Tn <- dim(theta_draws)[2]
Sn <- dim(theta_draws)[3]

# Sum by team-season
theta_med <- apply(theta_draws, c(2,3), median)
theta_lo  <- apply(theta_draws, c(2,3), quantile, probs = 0.025)
theta_hi  <- apply(theta_draws, c(2,3), quantile, probs = 0.975)

# labels 
season_levels <- sort(unique(df$season_id))  
team_levels   <- sort(unique(c(df$team_id_home, df$team_id_away)))

abbr_map <- unique(rbind(
  df[, c("team_id_home", "team_abbreviation_home")],
  setNames(df[, c("team_id_away", "team_abbreviation_away")],
           c("team_id_home","team_abbreviation_home"))
))
names(abbr_map) <- c("team_id","abbr")
abbr_map <- abbr_map[!duplicated(abbr_map$team_id), ]

# long data frame
plot_df <- expand_grid(team = 1:Tn, season = 1:Sn) %>%
  mutate(
    theta = as.vector(theta_med),
    lo    = as.vector(theta_lo),
    hi    = as.vector(theta_hi),
    season_id = factor(season_levels[season], levels = season_levels),
    team_id   = team_levels[team]
  ) %>%
  left_join(abbr_map, by = "team_id")


#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
# Alternative 2: Strength for season 
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################

logistic_2 <- "

  data {
  int<lower=1> G;
  int<lower=1> T;
  int<lower=1> S;

  array[G] int<lower=0,upper=1> y;
  array[G] int<lower=1,upper=T> home_team;
  array[G] int<lower=1,upper=T> away_team;
  array[G] int<lower=1,upper=S> season;
}

parameters {
  // season intercepts
  real alpha0;
  vector[S] alpha_raw;
  real<lower=0> sigma_alpha;

  // baseline team strength
  vector[T] theta0_raw;
  real<lower=0> sigma_theta0;

  // season deviations around baseline
  matrix[T,S] gamma_raw;
  real<lower=0> sigma_gamma;
}

transformed parameters {
  vector[S] alpha;
  vector[T] theta0;
  matrix[T,S] theta;   // theta_{t,s} = theta0_t + gamma_{t,s}

  // alpha_s centered
  alpha = alpha0 + sigma_alpha * (alpha_raw - mean(alpha_raw));

  // theta0 centered across teams
  theta0 = sigma_theta0 * (theta0_raw - mean(theta0_raw));

  // gamma centered within each season
  for (s in 1:S) {
    vector[T] g_s = to_vector(gamma_raw[, s]);
    theta[, s] = theta0 + sigma_gamma * (g_s - mean(g_s));
  }
}

model {
  // priors
  alpha0 ~ normal(0, 2.5);
  sigma_alpha  ~ normal(0, 1);
  sigma_theta0 ~ normal(0, 1);
  sigma_gamma  ~ normal(0, 1);

  alpha_raw  ~ normal(0, 1);
  theta0_raw ~ normal(0, 1);
  to_vector(gamma_raw) ~ normal(0, 1);

  // likelihood
  for (g in 1:G) {
    real eta =
      alpha[season[g]]
      + theta[home_team[g], season[g]]
      - theta[away_team[g], season[g]];

    y[g] ~ bernoulli_logit(eta);
  }
}

generated quantities {
  array[G] int y_rep;
  vector[G] log_lik;

  for (g in 1:G) {
    real eta =
      alpha[season[g]]
      + theta[home_team[g], season[g]]
      - theta[away_team[g], season[g]];

    y_rep[g]  = bernoulli_logit_rng(eta);
    log_lik[g]= bernoulli_logit_lpmf(y[g] | eta);
  }
}
"



sm_2 <- stan_model(model_code = logistic_2)
#saveRDS(sm_2, 'compiled_sm_2.rds')

fit_model2 <- sampling(
  object = sm_2,
  data = stan_data,
  seed = 9876,
  chains = 4,
  iter = 5000,     
  warmup = 1000,
  control = list(adapt_delta = 0.9, max_treedepth = 12)
)



post2 <- rstan::extract(fit_model2)

# 1) chequeo rápido: ¿existe theta?
stopifnot(!is.null(post$theta))

theta_arr <- post2$theta              # [iter, T, S]
I <- dim(theta_arr)[1]
T <- dim(theta_arr)[2]
S <- dim(theta_arr)[3]

season_labels <- if (exists("season_levels")) season_levels else paste0("S", 1:S)
team_labels   <- paste0("Team ", 1:T)  # o usa tus abbr si las tienes

# 2) resúmenes por (team, season)
theta_med <- apply(theta_arr, c(2,3), median, na.rm = TRUE)
theta_lo  <- apply(theta_arr, c(2,3), \(x) quantile(x, 0.025, na.rm=TRUE, names=FALSE))
theta_hi  <- apply(theta_arr, c(2,3), \(x) quantile(x, 0.975, na.rm=TRUE, names=FALSE))

# 3) a formato largo (sin líos de orden)
df_med <- as.data.frame(as.table(theta_med)) %>% rename(team = Var1, season = Var2, p50 = Freq)
df_lo  <- as.data.frame(as.table(theta_lo))  %>% rename(team = Var1, season = Var2, lo  = Freq)
df_hi  <- as.data.frame(as.table(theta_hi))  %>% rename(team = Var1, season = Var2, hi  = Freq)

theta_sum <- df_med %>%
  left_join(df_lo, by=c("team","season")) %>%
  left_join(df_hi, by=c("team","season")) %>%
  mutate(
    team   = factor(as.integer(team), levels = 1:T, labels = team_labels),
    season = factor(as.integer(season), levels = 1:S, labels = season_labels)
  )

# 4) debug mínimo (si aquí sale NA/0 filas, ya sabes por qué no pinta)
print(dim(theta_sum))
print(summary(theta_sum$p50))
stopifnot(nrow(theta_sum) > 0)
stopifnot(any(is.finite(theta_sum$p50)))

# 5) plot (líneas + CrI)
ggplot(theta_sum, aes(x = season, y = p50, group = team)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, group = team), alpha = 0.06) +
  geom_line(alpha = 0.35, linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Season", y = expression(theta[t,s]),
       title = "Team strength by season (median + 95% CrI)") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 16),
        axis.text  = element_text(size = 12),
        plot.title = element_text(size = 16))


###########################################
## Comparisons
###########################################

# Get log_lik 
ll0 <- loo::extract_log_lik(fit_model0, parameter_name = "log_lik", merge_chains = FALSE)
ll1 <- loo::extract_log_lik(fit_model1, parameter_name = "log_lik", merge_chains = FALSE)
ll2 <- loo::extract_log_lik(fit_model2, parameter_name = "log_lik", merge_chains = FALSE)

# r_eff mejora la precisión del LOO
#r1 <- loo::relative_eff(exp(ll1))
#r2 <- loo::relative_eff(exp(ll2))

#loo1 <- loo::loo(ll1, r_eff = r1)
#loo2 <- loo::loo(ll2, r_eff = r2)

#print(loo1)
#print(loo2)

# Comparación directa 
#loo::loo_compare(loo1, loo2)


waic0 <- loo::waic(ll0)
#waic1 <- loo::waic(ll1)
waic2 <- loo::waic(ll2)
#loo::loo_compare(waic1, waic2)
loo::loo_compare(waic0, waic2)
