library(agridat)
library(gnm)
library(betareg)
library(ggplot2)
library(viridis)
library(latex2exp)
library(rstan)
library(rstanarm)
library(coda)

# log-quasi-likelihood function for V(\mu) = \mu^d(1-\mu)^d
quasi_likelihood <- function(y, X, beta, psi, d)
{
  mu <- plogis(X %*% beta)
  require(cubature)
  y_mu <- cbind(y,mu)
  out <- apply(y_mu, 1, function(x) cubintegrate(f = quasi_score,
                                                 lower = x[1],
                                                 upper = x[2],
                                                 y = x[1],
                                                 d = d)$integral)
  sum(out) / psi
}

# gradient of log-quasi-likelihood with respct to \mu
quasi_score <- function(mu, y, d) (y - mu) / (mu^d * (1-mu)^d)

# gradient of log-quasi-likelihood with respect to \beta and V(\mu) = \mu^d(1-\mu)^d
gradient_quasi_likelihood <- function(y, X, beta, psi, d)
{
  mu <- plogis(X %*% beta)
  p <- length(beta)
  out <- (y - mu) / (psi * mu^(d-1) * (1-mu)^(d-1))
  out <- matrix(rep(out, p), ncol = p) * X
  colSums(out)
}

# log-quasi-posterior with V(\mu) = \mu^d(1-\mu)^d
logquasipost <- function(y, X, beta, psi, d)
{
  quasi_likelihood(y, X, beta, psi, d) + sum(dnorm(beta, 0, 10, log = TRUE))
}

# log-quasi-likelihood by Wedderburn in closed form
quasi_likelihood_w <- function(y, X, beta, psi)
{
  eta <- X %*% beta
  # mu <- exp(eta) / (1 + exp(eta))
  mu <- plogis(eta)
  sum((2*y + 1)*(log(mu)-log(1-mu)) - y/mu - (1-y)/(1-mu)) / psi
}

# Wedderburn log-quasi-posterior
logquasipost_w <- function(y, X, beta, psi)
{
  quasi_likelihood_w(y, X, beta, psi) + sum(dnorm(beta, 0, 10, log = TRUE))
}

# information matrix for Wedderburn quasi-posterior
information_ql_w <- function(y, X, beta, psi)
{
  eta <- X %*% beta
  t(X) %*% X / psi
}

# Random Walk Metropolis-Hastings algorithm for Wedderburn quasi-posterior
RWMH_ql_w <- function(R, burn_in, y, X, psi, start_beta, cov_init, eps, seed = 1)
{
  set.seed(seed)
  p <- length(start_beta)
  out <- matrix(0, nrow = R, ncol = p)
  beta <- start_beta
  logqp <- logquasipost_w(y, X, beta, psi)
  
  require(mvtnorm)
  require(matrixcalc)
  require(corpcor)
  sigma_rw <- eps*cov_init
  sigma_rw <- round(sigma_rw, 8)
  sigma_rw <- corpcor::make.positive.definite(sigma_rw, tol = 1e-8)
  
  for (r in 1:(burn_in + R))
  {
    beta_new <- c(beta + rmvnorm(1, sigma = sigma_rw))
    logqp_new <- logquasipost_w(y, X, beta_new, psi)
    alpha <- min(1, exp(logqp_new - logqp))
    if (runif(1) < alpha)
    {
      beta <- beta_new
      logqp <- logqp_new
    }
    if (r > burn_in)
    {
      out[r-burn_in,] <- beta
    }
    
    if (r %% 100 == 0) print(r)
  }
  
  out
}

# Random Walk Metropolis-Hastings algorithm for quasi-posterior with V(\mu) = \mu^d(1-\mu)^d
RWMH_ql <- function(R, burn_in, y, X, psi, d, start_beta, cov_init, eps, seed = 1)
{
  set.seed(seed)
  p <- length(start_beta)
  out <- matrix(0, nrow = R, ncol = p)
  beta <- start_beta
  logqp <- logquasipost(y, X, beta, psi, d)
  
  require(mvtnorm)
  require(matrixcalc)
  require(corpcor)
  sigma_rw <- eps*cov_init
  sigma_rw <- round(sigma_rw, 8)
  sigma_rw <- corpcor::make.positive.definite(sigma_rw, tol = 1e-8)
  
  for (r in 1:(burn_in + R))
  {
    beta_new <- c(beta + rmvnorm(1, sigma = sigma_rw))
    logqp_new <- logquasipost(y, X, beta_new, psi, d)
    alpha <- min(1, exp(logqp_new - logqp))
    if (runif(1) < alpha)
    {
      beta <- beta_new
      logqp <- logqp_new
    }
    if (r > burn_in)
    {
      out[r-burn_in,] <- beta
    }
    
    if (r %% 100 == 0) print(r)
  }
  
  out
}


# Data
library(agridat)
wed <- wedderburn.barley 
wed$y2 <- wed$y / 100 + 1e-7
y <- wed$y2
X <- model.matrix(~ gen + site, data = wed)
               
# Data as in Wedderburn 1974
data_book <- rbind(wed$y2[wed$site == "S1"], wed$y2[wed$site == "S2"], 
                   wed$y2[wed$site == "S3"], wed$y2[wed$site == "S4"],
                   wed$y2[wed$site == "S5"], wed$y2[wed$site == "S6"], 
                   wed$y2[wed$site == "S7"], wed$y2[wed$site == "S8"],
                   wed$y2[wed$site == "S9"])
mean_site <- rowMeans(data_book)
mean_gen <- colMeans(data_book)
var_site <- apply(data_book, 1, var)
var_gen <- apply(data_book, 2, var)

               
################################################################################
#           EXPLORATORY ANALYSIS FOR CHOICE OF VARIANCE FUNCTION               #
################################################################################
               
# Quasi-binomial setting
m_qb <- glm(y2 ~ gen + site, data = wed, family = quasibinomial)
p_qb <- predict(m_qb, type = "response")
pred_qb <- rbind(p_qb[wed$site == "S1"], p_qb[wed$site == "S2"], 
                 p_qb[wed$site == "S3"], p_qb[wed$site == "S4"],
                 p_qb[wed$site == "S5"], p_qb[wed$site == "S6"], 
                 p_qb[wed$site == "S7"], p_qb[wed$site == "S8"],
                 p_qb[wed$site == "S9"])
mean_site_qb <- rowMeans(pred_qb)
mean_gen_qb <- colMeans(pred_qb)
var_site_qb <- apply(pred_qb, 1, var)
var_gen_qb <- apply(pred_qb, 2, var)

# Beta regression setting
m_beta <- betareg(y2 ~ gen + site, data = wed, link = "logit", link.phi = "log")
p_beta <- predict(m_beta, type = "response")
pred_beta <- rbind(p_beta[wed$site == "S1"], p_beta[wed$site == "S2"], 
                   p_beta[wed$site == "S3"], p_beta[wed$site == "S4"],
                   p_beta[wed$site == "S5"], p_beta[wed$site == "S6"], 
                   p_beta[wed$site == "S7"], p_beta[wed$site == "S8"],
                   p_beta[wed$site == "S9"])
mean_site_beta <- rowMeans(pred_beta)
mean_gen_beta <- colMeans(pred_beta)
var_site_beta <- apply(pred_beta, 1, var)
var_gen_beta <- apply(pred_beta, 2, var)

# Wedderburn setting
m_wed <- glm(y2 ~ gen + site, data = wed, family = "wedderburn")
p_wed <- predict(m_wed, type = "response")
pred_wed <- rbind(p_wed[wed$site == "S1"], p_wed[wed$site == "S2"], 
                  p_wed[wed$site == "S3"], p_wed[wed$site == "S4"],
                  p_wed[wed$site == "S5"], p_wed[wed$site == "S6"], 
                  p_wed[wed$site == "S7"], p_wed[wed$site == "S8"],
                  p_wed[wed$site == "S9"])
mean_site_wed <- rowMeans(pred_wed)
mean_gen_wed <- colMeans(pred_wed)
var_site_wed <- apply(pred_wed, 1, var)
var_gen_wed <- apply(pred_wed, 2, var)

# Ad-hoc quasi-posterior setting with V(\mu) = \mu^{9/4}(1-\mu)^{9/4}
mle_wed <- summary(m_wed)$coefficients[,1]
loglikeq <- nlminb(start = mle_wed, objective = function(x) 
                      -quasi_likelihood(y = y, X = X, beta = x, psi = 1, d = 9/4))
beta_mle <- loglikeq$par
p_ql <- plogis(X %*% beta_mle)
pred_ql <- rbind(p_ql[wed$site == "S1"], p_ql[wed$site == "S2"], 
                 p_ql[wed$site == "S3"], p_ql[wed$site == "S4"],
                 p_ql[wed$site == "S5"], p_ql[wed$site == "S6"], 
                 p_ql[wed$site == "S7"], p_ql[wed$site == "S8"],
                 p_ql[wed$site == "S9"])
mean_site_ql <- rowMeans(pred_ql)
mean_gen_ql <- colMeans(pred_ql)
var_site_ql <- apply(pred_ql, 1, var)
var_gen_ql <- apply(pred_ql, 2, var)

# Comparison between empirical risk minimizers
v1 <- as.factor(c(as.character(wed$site), as.character(wed$gen)))
v2 <- rep(wed$y2, 2)
v3 <- as.factor(c(rep("Site", 90), rep("Variety", 90)))
v4 <- c(rep(mean_site_ql, each = 10), rep(mean_gen_ql, 9))
v5 <- c(rep(mean_site_qb, each = 10), rep(mean_gen_qb, 9))
v6 <- c(rep(mean_site_wed, each = 10), rep(mean_gen_wed, 9))
v7 <- c(rep(mean_site_beta, each = 10), rep(mean_gen_beta, 9))
wed_df <- data.frame(v1, v2, v3, v4, v5, v6, v7)
mean_plot <- ggplot(wed_df) +
  geom_point(aes(x = v1, y = v2, shape = "obs"), size = 2) + 
  facet_grid(~ v3, scales = "free_x") +
  geom_point(aes(x = v1, y = v4, col = "ql", shape = "obs"), size = 2) +
  geom_point(aes(x = v1, y = v5, col = "bin", shape = "obs"), size = 2) +
  geom_point(aes(x = v1, y = v6, col = "wed", shape = "obs"), size = 2) +
  geom_point(aes(x = v1, y = v7, col = "beta", shape = "obs"), size = 2) +
  stat_summary(aes(x = v1, y = v4, group = v3, col = "ql"), 
               fun = "mean", geom = "line") +
  stat_summary(aes(x = v1, y = v5, group = v3, col = "bin"), 
               fun = "mean", geom = "line") +
  stat_summary(aes(x = v1, y = v6, group = v3, col = "wed"), 
               fun = "mean", geom = "line") +
  stat_summary(aes(x = v1, y = v7, group = v3, col = "beta"), 
               fun = "mean", geom = "line") +
  scale_color_manual(name = "Assumed variance",
                     values = c("beta" = "#0033CC",
                                "bin" = "orange",
                                "wed" = "#339933",
                                "ql" = "#CC3300"),
                     labels = c(TeX(r'($var(Y) = \mu(1-\mu)(\phi+1)^{-1}$)'),
                                TeX(r'($var(Y) = \psi\mu(1-\mu)$)'),
                                TeX(r'($var(Y) = \psi\mu^2(1-\mu)^2$)'),
                                TeX(r'($var(Y) = \psi\mu^{9/4}(1-\mu)^{9/4}$)')
                     )) +
  scale_shape_manual(name = "",
                     values = c("obs" = 19), 
                     labels = c("Observed values")) +
  labs(x = "", y = "Y") +
  theme_light() +
  theme(legend.text.align = 0,
        strip.text = element_text(size = 14, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"),
        legend.position = "top",
        legend.text = element_text(size = 11))
mean_plot

# Conparison theoretical variance function vs observed
mu_val <- seq(0.00001, 0.99999, length = 1000)
v <- function(mu, p, psi = 1) psi * mu^p * (1-mu)^p
v_beta <- function(mu, v) mu * (1-mu) / (v + 1)
mu_i <- plogis(X %*% beta_mle)
p <- length(beta_mle)
N <- length(y)
d <- 9/4
psi <- sum( (y - mu_i)^2 / (mu_i^d * (1-mu_i)^d) ) / (N - p)
df_ql <- data.frame(x = mu_val, y = v(mu_val, p = d, psi = psi))
psi_wed <- summary(m_wed)$dispersion
df_wed <- data.frame(x = mu_val, y = v(mu_val, p = 2, psi = psi_wed))
psi_bin <- summary(m_qb)$dispersion
df_bin <- data.frame(x = mu_val, y = v(mu_val, p = 1, psi = psi_bin))
prec_beta <- summary(m_beta)$coefficients$precision[1]
df_beta <- data.frame(x = mu_val, y = v_beta(mu_val, v = prec_beta))
df_var <- data.frame(mu = c(mean_site, mean_gen),
                     v = c(var_site, var_gen),
                     group = c(rep("Site", length(mean_site)), 
                               rep("Variety", length(mean_gen))))
df_var$group <- as.factor(df_var$group)
var_plot <- ggplot(data = df_var, aes(x = mu, y = v)) +
  geom_point(aes(shape = "obs"), size = 2.7, color = "black", fill = "black") +
  facet_wrap(~ group) +
  geom_line(data = df_ql, aes(x = x, y = y, color = "ql")) +
  geom_line(data = df_wed, aes(x = x, y = y, color = "wed")) +
  geom_line(data = df_bin, aes(x = x, y = y, color = "bin")) +
  geom_line(data = df_beta, aes(x = x, y = y, color = "beta")) +
  scale_color_manual(name = "Assumed variance",
                     values = c("beta" = "#0033CC",
                                "bin" = "orange",
                                "wed" = "#339933",
                                "ql" = "#CC3300"),
                     labels = c(TeX(r'($var(Y) = \mu(1-\mu)(\phi+1)^{-1}$)'),
                                TeX(r'($var(Y) = \psi\mu(1-\mu)$)'),
                                TeX(r'($var(Y) = \psi\mu^2(1-\mu)^2$)'),
                                TeX(r'($var(Y) = \psi\mu^{9/4}(1-\mu)^{9/4}$)')
                     )) +
  scale_shape_manual(name = "",
                     values = c("obs" = 19), 
                     labels = c("Observed values")) +
  labs(x = "Mean", y = "Variance") +
  theme_light() + 
  theme(legend.text.align = 0,
        strip.text = element_text(size = 14, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"),
        legend.position = "top",
        legend.text = element_text(size = 11))
var_plot


################################################################################
#                                   MODELING                                   #
################################################################################
                   
# BETA Regression
GLM_stan_BB <- stan_betareg(y ~ X - 1,
                            prior = normal(0, 10),
                            link = "logit",
                            link.phi = "log",
                            chains = 4,
                            cores = 4,
                            iter = 5000,
                            seed = 15)
GLM_out_BB <- as.mcmc(as.matrix(GLM_stan_BB))
plot(GLM_out_BB)
beta_bb <- colMeans(GLM_out_BB[,-NCOL(GLM_out_BB)])
beta_phi <- mean(GLM_out_BB[,NCOL(GLM_out_BB)])

# QUASI-LIKELIHOOD GLM d = 9/4
S <- solve(optimHess(beta_mle, function(x) -logquasipost(y = y, X = X, beta = x,
                                                         psi = psi, d = d)))
GLM_out <- as.mcmc(RWMH_ql(R = 70000, burn_in = 5000,
                           y = y, X = X, psi = psi, d = d,
                           start_beta = beta_mle, cov_init = S, 
                           eps = (2.38^2)/p, seed = 1))
keep <- seq(1:10000) * 7 
GLM_thin <- as.mcmc(GLM_out[keep,])
plot(GLM_thin)
beta_ql <- colMeans(GLM_thin)

# QUASI-LIKELIHOOD GLM d = 2 (Wedderburn)
beta_mle_w <- as.numeric(summary(m_wed)$coef[,1])
S_w <- solve(information_ql_w(y, X, beta_mle_w, psi_wed))
GLM_out_w <- as.mcmc(RWMH_ql_w(R = 70000, burn_in = 5000,
                               y = y, X = X, psi = psi_wed,
                               start_beta = beta_mle_w, cov_init = S_w,
                               eps = (2.38^2)/p, seed = 2))
keep_w <- seq(1:10000) * 7
GLM_thin_w <- as.mcmc(GLM_out_w[keep_w,])
plot(GLM_thin_w)
beta_ql_w <- colMeans(GLM_thin_w)


################################################################################
#                               GOODNESS-OF-FIT                                #
################################################################################
                     
idx_1 <- 42; site_1 <- max(X[idx_1, c(1,11:18)] * c(1:9)); gen_1 <- max(X[idx_1, 1:10] * c(1:10))
idx_2 <- 70; site_2 <- max(X[idx_2, c(1,11:18)] * c(1:9)); gen_2 <- max(X[idx_2, 1:10] * c(1:10))
idx_3 <- 88; site_3 <- max(X[idx_3, c(1,11:18)] * c(1:9)); gen_3 <- max(X[idx_3, 1:10] * c(1:10))
df_prev <- data.frame(x_ql = c(plogis(GLM_thin %*% X[idx_1,]),
                               plogis(GLM_thin %*% X[idx_2,]),
                               plogis(GLM_thin %*% X[idx_3,])),
                      x_bb = c(plogis(GLM_out_BB[,-NCOL(GLM_out_BB)] %*% X[idx_1,]),
                               plogis(GLM_out_BB[,-NCOL(GLM_out_BB)] %*% X[idx_2,]),
                               plogis(GLM_out_BB[,-NCOL(GLM_out_BB)] %*% X[idx_3,])),
                      x_w = c(plogis(GLM_thin_w %*% X[idx_1,]),
                              plogis(GLM_thin_w %*% X[idx_2,]),
                              plogis(GLM_thin_w %*% X[idx_3,])),
                      obs = factor(c(rep(paste0("Site ",site_1,", Variety ",gen_1), 1e4),
                                     rep(paste0("Site ",site_2,", Variety ",gen_2), 1e4),
                                     rep(paste0("Site ",site_3,", Variety ",gen_3), 1e4))),
                      v_line = c(rep(data_book[site_1, gen_1], 1e4),
                                 rep(data_book[site_2, gen_2], 1e4),
                                 rep(data_book[site_3, gen_3], 1e4)))
prev_plot <- ggplot(df_prev) +
  facet_wrap(~ obs, scales = "free") +
  stat_density(aes(x = x_ql, col = "ql"),
               geom = "line", position = "identity", adjust = 2) +
  stat_density(aes(x = x_w, col = "wed"),
               geom = "line", position = "identity", adjust = 2) +
  stat_density(aes(x = x_bb, col = "beta"),
               geom = "line", position = "identity", adjust = 2) +
  scale_color_manual(name = "Assumed variance",
                     values = c("beta" = "#0033CC",
                                "wed" = "#339933",
                                "ql" = "#CC3300"),
                     labels = c(TeX(r'($var(Y) = \mu(1-\mu)(\phi+1)^{-1}$)'),
                                TeX(r'($var(Y) = \psi\mu^2(1-\mu)^2$)'),
                                TeX(r'($var(Y) = \psi\mu^{9/4}(1-\mu)^{9/4}$)')
                     )) +
  geom_vline(aes(xintercept = v_line, linetype = "obs")) +
  scale_linetype_manual(name = "",
                        values = c("obs" = "dashed"),
                        labels = c("Observed value")) +
  labs(x = TeX(r'($\mu$)'), y = "Density") +
  theme_light() + 
  guides(color = guide_legend(order = 1),
         linetype = guide_legend(order = 2)) +
  theme(legend.text.align = 0,
        strip.text = element_text(size = 12, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"),
        legend.position = "top",
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        legend.text = element_text(size = 11))

