library(rstan)
library(rstanarm)
library(coda)
library(ggplot2)
library(latex2exp)

n_logqpost <- function(beta, y, X) {
  mu <- X %*% beta
  - sum((mu - y + 1)*exp(-mu))
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

beta_true <- c(-3, 2, 1.5, 1)
psi_true <- 2.5
N <- 300
p <- length(beta_true)
S <- 100  
seed_sim <- 2023
beta_mle <- matrix(NA, nrow = S, ncol = p)
mu_ql <- mu_lm <- matrix(NA, nrow = S, ncol = N)
psi <- rep(NA, S)
y_sim <- X_sim <- list()
beta_ql <- beta_lm <- sigma <- freq_est <- list()

for (s in 1:S) {
  
  print(s)
  
  # Simulated data
  set.seed(seed_sim + s)
  X <- cbind(rep(1, N), rnorm(N, mean = 0, sd = 1), 
             rnorm(N, mean = 0, sd = 1), rnorm(N, mean = 0, sd = 1))
  X <- as.matrix(X)
  mu_true = X %*% beta_true
  var_true = psi_true * exp(mu_true)
  y <- rnorm(N, mean = mu_true, sd = sqrt(var_true))
  X_sim[[s]] <- X
  y_sim[[s]] <- y
  
  # Linear model
  LM_stan <- stan_glm(y ~ X - 1,
                      family = gaussian,
                      prior = normal(0, 10),
                      chains = 3,
                      cores = 3,
                      warmup = 1000,
                      iter = 3000)
  beta_lm[[s]] <- as.mcmc(as.matrix(LM_stan)[,1:p])
  sigma[[s]] <- as.mcmc(as.matrix(LM_stan)[,p+1])
  
  # Quasi-posterior
  beta_mle[s,] <- nlminb(start = rep(0, p), 
                         objective = n_logqpost,
                         y = y, X = X)$par
  fit <- X %*% beta_mle[s,]
  psi[s] <- sum( (y - fit)^2 / exp(fit) ) / (N - p)
  stan_data <- as.list(as.data.frame(X))
  names(stan_data) <- c("intercept", "v1", "v2", "v3")
  stan_data$N <- N
  stan_data$p <- p
  stan_data$y <- y
  stan_data$rec_psi <- 1 / psi[s]
  fit <- stan(file = "sim_heteroscedastic_data.stan",
              data = stan_data,
              warmup = 1000,
              iter = 3000, 
              chains = 3, 
              cores = 3)
  GLM_ql <- as.mcmc(as.matrix(fit))
  beta_ql[[s]] <- GLM_ql[,1:stan_data$p]
  
}
# Important: check the convergence of MCMC

# Coverage
cover <- cbind(lm_90 = coverage(beta_lm, beta_true, 0.90),
               ql_90 = coverage(beta_ql, beta_true, 0.90),
               lm_95 = coverage(beta_lm, beta_true, 0.95),
               ql_95 = coverage(beta_ql, beta_true, 0.95),
               lm_99 = coverage(beta_lm, beta_true, 0.99),
               ql_99 = coverage(beta_ql, beta_true, 0.99))
rownames(cover) <- c("beta_1", "beta_2", "beta_3", "beta_4")
cover

# Figure
df <- data.frame(x = factor(rep(c("Quasi-posterior", "Linear model"), each = p*S)),
                 y = c(t(sapply(beta_ql, colMeans)), t(sapply(beta_lm, colMeans))),
                 beta = rep(rep(1:p, each = S), 2))
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
