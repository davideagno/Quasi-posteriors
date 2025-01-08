rm(list = ls())

library(rstan)
library(rstanarm)
library(coda)
library(ggplot2)
library(latex2exp)
library(TeachingDemos)

coverage <- function(beta, beta_true, level) {
  emp_hpd_lower <- function(x) TeachingDemos::emp.hpd(x, conf = level)[1]
  emp_hpd_upper <- function(x) TeachingDemos::emp.hpd(x, conf = level)[2]
  
  S <- length(beta); p <- length(beta_true)
  hpd_low <- hpd_upp <- matrix(NA, nrow = S, ncol = p)
  
  for(s in 1:S) {
    hpd_low[s,] <- apply(beta[[s]], 2, emp_hpd_lower)
    hpd_upp[s,] <- apply(beta[[s]], 2, emp_hpd_upper)
  }
  
  out <- rep(NA, p)
  for (r in 1:p) {
    out[r] <- sum((hpd_low[,r] < beta_true[r]) & (beta_true[r] < hpd_upp[,r])) / S
  }
  
  out
}

boot_psi <- function(y, X, y_B, X_B) {
  mod <- glm(y_B ~ X_B - 1, family = quasipoisson)
  beta_mle_B <- mod$coefficients
  pred <- exp(c(X %*% beta_mle_B))
  dbl <- colSums((y - pred) * X)
  dbp <- beta_mle_B / 6.25
  d2 <- t(X) %*% diag(pred) %*% X
  num <- as.numeric(dbl %*% dbp) + sum(diag(d2))
  den <- sum(dbl^2)
  return(list(num = num, den = den))
}
psi_bmk_par <- function(y, X, B, seed = 123, cores = 7) {
  require(parallel)
  nn <- length(y)
  idx_boot <- matrix(nrow = nn, ncol = B)
  for (b in 1:B) {
    set.seed(seed + b)
    idx_boot[,b] <- sample(1:nn, size = nn, replace = TRUE)
  }
  
  mclapply_function <- function(i) {
    boot_psi(y, X, y[idx_boot[,i]], X[idx_boot[,i],])
  }
  out_boot_psi <- mclapply(1:B, mclapply_function, mc.cores = cores)
  
  out <- matrix(unlist(out_boot_psi), ncol = B)
  sum(out[1,]) / sum(out[2,])
}
psi_bmk <- function(y, X, B, seed = 123) {
  num <- den <- rep(NA, B)
  nn <- length(y)
  for (b in 1:B) {
    set.seed(seed + b)
    idx_boot <- sample(1:nn, size = nn, replace = TRUE)
    y_boot <- y[idx_boot]
    X_boot <- X[idx_boot,]
    bb <- boot_psi(y, X, y_boot, X_boot)
    num[b] <- bb$num
    den[b] <- bb$den
  }
  sum(num) / sum(den)
}

dfd <- function(beta, y, X) {
  mu <- exp(X %*% beta)
  sum((y^2) / (mu^2) - 2*(y + 1)/mu)
}
boot_psi_dfd <- function(y, X, y_B, X_B) {
  mm <- nlminb(start = rep(1, NCOL(X)), objective = dfd, y = y_B, X = X_B)
  beta_hat_B <- c(mm$par)
  pred <- exp(c(X %*% beta_hat_B))
  pred[pred < 1e-150] <- 1e-150
  dbl <- colSums((-2*(y^2)/(pred^2) + 2*(y+1)/pred) * X)
  dbp <- beta_hat_B / 6.25
  d2 <- t(X) %*% diag(4*(y^2)/(pred^2) - 2*(y+1)/pred) %*% X
  num <- as.numeric(dbl %*% dbp) + sum(diag(d2))
  den <- sum(dbl^2)
  return(list(num = num, den = den))
}
psi_dfd_par <- function(y, X, B, seed = 123, cores = 7) {
  require(parallel)
  nn <- length(y)
  idx_boot <- matrix(nrow = nn, ncol = B)
  for (b in 1:B) {
    set.seed(seed + b)
    idx_boot[,b] <- sample(1:nn, size = nn, replace = TRUE)
  }
  
  mclapply_function <- function(i) {
    boot_psi_dfd(y, X, y[idx_boot[,i]], X[idx_boot[,i],])
  }
  out_boot_psi <- mclapply(1:B, mclapply_function, mc.cores = cores)
  
  out <- matrix(unlist(out_boot_psi), ncol = B)
  idx_out <- which(is.finite(out[2,]) & out[2,] < 1e6)
  sum(out[1, idx_out], na.rm = TRUE) / sum(out[2, idx_out], na.rm = TRUE)
}

# SETTINGS
beta_true <- c(3.5, 1.5, -1, 0.5)
psi_true <- 3.5
# N <- 1000
N <- 100
p <- length(beta_true)
S <- 200
seed_sim <- 2023
beta_mle <- matrix(NA, nrow = S, ncol = p)
mu_ql <- mu_poi <- mu_nb <- matrix(NA, nrow = S, ncol = N)
psi <- psi_bb <- ls_dfd <- rep(NA, S)
y_sim <- X_sim <- list()
beta_ql <- beta_poi <- beta_nb <- phi_nb <- freq_est <- beta_ql_bb <- beta_dfd <- list()
err_mean <- err_var <- matrix(NA, nrow = S, ncol = N)

# SIMULATE DATA
for(s in 1:S) {
  set.seed(seed_sim + s)
  X <- cbind(rep(1, N), rnorm(N, mean = 0, sd = 1),
             rnorm(N, mean = 0, sd = 1), rnorm(N, mean = 0, sd = 1))
  X <- as.matrix(X)
  mu_true <- exp(X %*% beta_true)
  var_true <- psi_true * mu_true
  y_star <- rgamma(N, shape = mu_true / psi_true, rate = 1 / psi_true)
  y <- round(y_star)
  X_sim[[s]] <- X
  y_sim[[s]] <- y
}

# QUASI POSTERIOR 
for (s in 1:S) {
  set.seed(s)
  print(s)
  X <- X_sim[[s]]
  y <- y_sim[[s]]
  mod <- glm(y ~ X - 1, family = quasipoisson)
  freq_est[[s]] <- summary(mod)$coefficients[,1:2]
  beta_mle[s,] <- mod$coefficients
  psi[s] <- summary(mod)$dispersion
  stan_data <- list()
  stan_data$X <- as.matrix(X)
  stan_data$N <- N
  stan_data$p <- p
  stan_data$y <- as.vector(y)
  stan_data$rec_psi <- 1 / psi[s]
  fit <- stan(file = "sim_counts_mat.stan", data = stan_data,
              warmup = 1000, iter = 2000,
              chains = 3, cores = 3)
  GLM_ql <- as.mcmc(as.matrix(fit))
  beta_ql[[s]] <- GLM_ql[,1:stan_data$p]
}

# QUASI POSTERIOR - BOOTSTRAP TUNING
for (s in 1:S) {
  set.seed(s)
  print(s)
  X <- X_sim[[s]]
  y <- y_sim[[s]]
  psi_bb[s] <- 1/psi_bmk_par(y, X, 500)
  stan_data <- list()
  stan_data$X <- as.matrix(X)
  stan_data$N <- N
  stan_data$p <- p
  stan_data$y <- as.vector(y)
  stan_data$rec_psi <- 1 / psi_bb[s]
  fit <- stan(file = "sim_counts_mat.stan", data = stan_data,
              warmup = 1000, iter = 2000,
              chains = 3, cores = 3)
  GLM_ql_bb <- as.mcmc(as.matrix(fit))
  beta_ql_bb[[s]] <- GLM_ql_bb[,1:stan_data$p]
}

# POISSON
for (s in 1:S) {
  set.seed(s)
  print(s)
  X <- X_sim[[s]]
  y <- y_sim[[s]]
  GLM_poi <- stan_glm(y ~ X - 1, family = poisson,
                      warmup = 1000, iter = 2000,
                      chains = 3, cores = 3)
  beta_poi[[s]] <- as.mcmc(as.matrix(GLM_poi))
}

# NEGATIVE BINOMIAL
for (s in 1:S) {
  print(s)
  X <- X_sim[[s]]
  y <- y_sim[[s]]
  GLM_nb <- stan_glm(y ~ X - 1, family = neg_binomial_2,
                     warmup = 1000, iter = 2000,
                     chains = 3, cores = 3)
  beta_nb[[s]] <- as.mcmc(as.matrix(GLM_nb)[,1:p])
}

# DFD-BAYES
for (s in 1:S) {
  set.seed(s)
  print(s)
  X <- X_sim[[s]]
  y <- y_sim[[s]]
  stan_data <- list()
  stan_data$X <- as.matrix(X)
  stan_data$N <- N
  stan_data$p <- p
  stan_data$y <- as.vector(y)
  ls_dfd[s] <- psi_dfd_par(y, X, 400)
  stan_data$rec_psi <- ls_dfd[s]
  fit <- stan(file = "sim_counts_dfd.stan", data = stan_data,
              warmup = 1000, iter = 2000, chains = 3, cores = 3)
  GLM_dfd <- as.mcmc(as.matrix(fit))
  beta_dfd[[s]] <- GLM_dfd[,1:stan_data$p]
  
}

post_mean <- cbind(poi = rowMeans(sapply(beta_poi, colMeans)),
                   nb = rowMeans(sapply(beta_nb, colMeans)),
                   dfd = rowMeans(sapply(beta_dfd, colMeans)),
                   ql = rowMeans(sapply(beta_ql, colMeans)),
                   ql_bt = rowMeans(sapply(beta_ql_bb, colMeans)))
rownames(post_mean) <- c("beta_1", "beta_2", "beta_3", "beta_4")
round(post_mean, 2)

sd_post_mean <- cbind(poi = apply(sapply(beta_poi, colMeans), 1, sd),
                      nb = apply(sapply(beta_nb, colMeans), 1, sd),
                      dfd = apply(sapply(beta_dfd, colMeans), 1, sd),
                      ql = apply(sapply(beta_ql, colMeans), 1, sd),
                      ql_bt = apply(sapply(beta_ql_bb, colMeans), 1, sd))
rownames(sd_post_mean) <- c("beta_1", "beta_2", "beta_3", "beta_4")
round(sd_post_mean, 3)

cover <- cbind(poi = coverage(beta_poi, beta_true, 0.95),
               nb = coverage(beta_nb, beta_true, 0.95),
               dfd = coverage(beta_dfd, beta_true, 0.95),
               ql = coverage(beta_ql, beta_true, 0.95),
               ql_bt = coverage(beta_ql_bb, beta_true, 0.95))
rownames(cover) <- c("beta_1", "beta_2", "beta_3", "beta_4")
round(cover, 3)


