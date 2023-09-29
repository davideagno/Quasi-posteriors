library(rstan)
library(rstanarm)
library(coda)
library(ggplot2)
library(latex2exp)
library(TeachingDemos)

bias_m1 <- function(shape, rate, r, a = 0) {
  sigma <- sqrt(shape) / rate
  serie <- rep(0, 200)
  den_const <- r*sqrt(shape)
  for (k in 1:length(serie)) {
    f1 <- (-1)^k / k
    frac <- 2*pi*k / den_const
    f2 <- (1 + frac^2)^(-shape/2)
    theta <- atan(frac)
    f3 <- sin(shape*theta - 2*pi*k*a)
    serie[k] <- f1*f2*f3
  }
  sum(serie) * r * sigma / pi
}

bias_m2 <- function(shape, rate, r, a = 0){
  sigma2 <- shape / rate^2
  serie <- rep(0, 200)
  den_const <- r*sqrt(shape)
  for (k in 1:length(serie)) {
    f1 <- (-1)^k
    frac <- 2*pi*k / den_const
    f2 <- (1 + frac^2)^(-shape/2)
    theta <- atan(frac)
    f3 <- cos(shape*theta - 2*pi*k*a) / (pi*k)
    f4 <- ((1 + frac^2) / shape)^(-0.5) * 2 / r
    f5 <- sin((shape + 1)*theta - 2*pi*k*a)
    serie[k] <- f1 * f2 * (f3 + f4) * f5
  }
  a1 <- r^2 / 12
  a2 <- r^2 / pi
  sigma2 * (a1 + a2*sum(serie))
}

ratio_var <- function(shape, rate, r, a = 0) {
  b2 <- bias_m2(shape, rate, r, a)
  bm <- bias_m1(shape, rate, r, a)
  M <- shape / rate
  V <- shape / rate^2
  bv <- b2 - bm^2 + 2*M*bm
  (V + bv) / V
}

coverage <- function(beta, beta_true, level)
{
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

beta_true <- c(3.5, 1.5, -1, 0.5)
psi_true <- 3.5
N <- 1000
p <- length(beta_true)
S <- 100 
seed_sim <- 2023
beta_mle <- matrix(NA, nrow = S, ncol = p)
mu_ql <- mu_poi <- mu_nb <- matrix(NA, nrow = S, ncol = N)
psi <- rep(NA, S)
y_sim <- X_sim <- list()
beta_ql <- beta_poi <- beta_nb <- phi_nb <- freq_est <- list()
err_mean <- err_var <- matrix(NA, nrow = S, ncol = N)

for (s in 1:S) {
  
  print(s)
  
  # Simulate data
  set.seed(seed_sim + s)
  X <- cbind(rep(1, N), rnorm(N, mean = 0, sd = 1),
             rnorm(N, mean = 0, sd = 1), rnorm(N, mean = 0, sd = 1))
  X <- as.matrix(X)
  mu_true <- exp(X %*% beta_true)
  var_true <- psi_true * mu_true
  for (i in 1:N) {
    err_var[s,i] <- ratio_var(shape = mu_true[i] / psi_true,
                              rate = 1 / psi_true,
                              r = 1 / sqrt(psi_true*mu_true[i]))
    err_mean[s,i] <- bias_m1(shape = mu_true[i] / psi_true,
                             rate = 1 / psi_true,
                             r = 1 / sqrt(psi_true*mu_true[i]))
  }
  y_star <- rgamma(N, shape = mu_true / psi_true, rate = 1 / psi_true)
  y <- round(y_star)
  X_sim[[s]] <- X
  y_sim[[s]] <- y
  
  # Poisson model
  GLM_poi <- stan_glm(y ~ X - 1, family = poisson,
                      warmup = 1000, iter = 2000,
                      chains = 3, cores = 3)
  beta_poi[[s]] <- as.mcmc(as.matrix(GLM_poi))
  
  # Negative-binomial model
  GLM_nb <- stan_glm(y ~ X - 1, family = neg_binomial_2,
                     warmup = 1000, iter = 2000,
                     chains = 3, cores = 3)
  beta_nb[[s]] <- as.mcmc(as.matrix(GLM_nb)[,1:p])
  
  # Quasi-posterior
  mod <- glm(y ~ X - 1, family = quasipoisson)
  freq_est[[s]] <- summary(mod)$coefficients[,1:2]
  beta_mle[s,] <- mod$coefficients
  psi[s] <- summary(mod)$dispersion
  stan_data <- as.list(as.data.frame(X))
  names(stan_data) <- c("intercept", "v1", "v2", "v3")
  stan_data$N <- N
  stan_data$p <- p
  stan_data$y <- y
  stan_data$rec_psi <- 1 / psi[s]
  fit <- stan(file = "sim_counts.stan", data = stan_data,
              warmup = 1000, iter = 2000,
              chains = 3, cores = 3)
  GLM_ql <- as.mcmc(as.matrix(fit))
  beta_ql[[s]] <- GLM_ql[,1:stan_data$p]
  
}
# Important: check the MCMC convergence

# Coverage
cover <- cbind(nb_90 = coverage(beta_nb, beta_true, 0.90),
               poi_90 = coverage(beta_poi, beta_true, 0.90),
               ql_90 = coverage(beta_ql, beta_true, 0.90),
               nb_95 = coverage(beta_nb, beta_true, 0.95),
               poi_95 = coverage(beta_poi, beta_true, 0.95),
               ql_95 = coverage(beta_ql, beta_true, 0.95),
               nb_99 = coverage(beta_nb, beta_true, 0.99),
               poi_99 = coverage(beta_poi, beta_true, 0.99),
               ql_99 = coverage(beta_ql, beta_true, 0.99))
rownames(cover) <- c("beta_1", "beta_2", "beta_3", "beta_4")
cover

# Figure
df <- data.frame(x = factor(rep(c("QP", "POI", "NB"), each = p*S)),
                 y = c(t(sapply(beta_ql, colMeans)), 
                       t(sapply(beta_poi, colMeans)), 
                       t(sapply(beta_nb, colMeans))),
                 beta = rep(rep(1:p, each = S), 3))
df_true <- data.frame(val = beta_true, beta = 1:p)
fig <- ggplot(df, aes(x = x, y = y)) +
  facet_wrap(~ beta, scales = "free", nrow = 1, 
             labeller = label_bquote(beta[.(beta)])) +
  geom_boxplot(fill = "gray85", outlier.size = 1) + 
  geom_hline(data = df_true, aes(yintercept = val, linetype = "true"), 
             color = "red") +
  scale_linetype_manual(name = "",
                        values = c("true" = "dashed"),
                        labels = c("True value")) +
  labs(x = "", y = "Posterior mean") +
  theme_light() + 
  theme(legend.position = "top",
        strip.text = element_text(size = 12, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"))
