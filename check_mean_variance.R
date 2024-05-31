rm(list = ls())
library(MASS)
library(ggplot2)
library(latex2exp)

beta_true <- c(2, 1)
psi_true <- 2.5
N <- 500
set.seed(123)
X <- as.matrix(cbind(rep(1, N), rnorm(N, 0, 1)))
mu_true <- exp(X %*% beta_true)
var_true <- psi_true * mu_true
shape_true <- mu_true^2 / var_true
rate_true <- mu_true / var_true
y_star <- rgamma(N, shape = shape_true, rate = rate_true)
y <- round(y_star)


################################################
####               CHECK MEAN               ####
################################################
lm_or <- lm(y ~ X - 1)$coefficients
log_y <- log(y)
lm_log <- lm(log_y[log_y > -Inf] ~ X[log_y > -Inf,] - 1)$coefficients
lm_sqrt <- lm(sqrt(y) ~ X - 1)$coefficients
lm_rcub <- lm(y^(1/3) ~ X - 1)$coefficients

log_y_plot <- log_y[log_y > -Inf]
df <- data.frame(x = c(rep(X[,2], 3), X[log_y > -Inf, 2]),
                 y = c(y, sqrt(y), y^(1/3), log_y_plot),
                 group = c(rep(c(1,3,4), each = length(y)), rep(c(2), length(log_y_plot))))
df$group <- as.factor(df$group)
df_lm <- data.frame(inte = c(lm_or[1], lm_sqrt[1], lm_rcub[1], lm_log[1]),
                    slop = c(lm_or[2], lm_sqrt[2], lm_rcub[2], lm_log[2]),
                    group = as.factor(c(1,3,4,2)))
levels(df$group) <- levels(df_lm$group) <- c("Identity", "Logarithm", "Square root", "Cubic root")
mean_plot <- ggplot(data = df, aes(x = x, y = y)) +
  geom_point(size = 1, color = "#336666") +
  facet_wrap(~ group, scales = "free_y", nrow = 1) +
  geom_abline(data = df_lm, aes(intercept = inte, slope = slop), color = "#CC3300") +
  labs(x = "X", y = "") +
  theme_light() +
  theme(strip.text = element_text(size = 12, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"))
mean_plot


####################################################
####               CHECK VARIANCE               ####
####################################################

# Poisson
poi <- glm(y ~ X - 1, family = poisson)
pred_poi <- fitted(poi)

# Negative binomial
nb <- glm.nb(y ~ X - 1)
(phi_nb <- nb$theta)
pred_nb <- fitted(nb)

# V(\mu) = \mu (quasi-Poisson)
qp <- glm(y ~ X - 1, family = quasipoisson)
(psi_qp <- summary(qp)$dispersion)
beta_qp <- c(summary(qp)$coefficients[,1])
pred_qp <- fitted(qp)

# V(\mu) = \mu^2 (quasi-Gamma)
nlog_qpost_mle_2 <- function(beta, y, X) {
  mu <- exp(X %*% beta)
  - sum( - y / mu - log(mu))
}
qp2_mle <- nlminb(c(0,0), function(b) nlog_qpost_mle_2(b, y, X))
qp2_mle
beta_qp2 <- qp2_mle$par
pred_qp2 <- c(exp(X %*% beta_qp2))
(psi_qp2 <- sum((y - pred_qp2)^2 / (pred_qp2^2)) / (N - 2))

# V(\mu) = exo(\mu)
nlog_qpost_mle_exp <- function(beta, y, X) {
  mu <- exp(X %*% beta)
  sum((-mu + y - 1)*exp(-mu))
}
qp_exp_mle <- nlminb(c(0,0), function(b) nlog_qpost_mle_exp(b, y, X))
qp_exp_mle
beta_qp_exp <- qp_exp_mle$par
pred_qp_exp <- c(exp(X %*% beta_qp_exp))
(psi_qp_exp <- sum((y - pred_qp_exp)^2 / (exp(pred_qp_exp))) / (N - 2))


x_class <- round(X[,2], 1)
range_val <- round(seq(-1.5, 1.5, by = 0.1), 1)
mean_obs <- var_obs <- var_true <- var_poi <- var_nb <- var_qp <- var_qp2 <- var_qp_exp <- rep(NA, length(range_val))
for (i in 1:length(range_val)) {
  var_true[i] <- psi_true * exp(sum(beta_true * c(1, range_val[i])))^2
  mean_obs[i] <- mean(y[x_class == range_val[i]])
  var_obs[i] <- var(y[x_class == range_val[i]])
  var_poi[i] <- mean_obs[i]
  var_nb[i] <- mean_obs[i]*(1 + mean_obs[i]/phi_nb)
  var_qp[i] <- psi_qp * mean_obs[i]
  var_qp2[i] <- psi_qp2 * mean_obs[i]^2
  var_qp_exp[i] <- psi_qp_exp * exp(mean_obs[i])
}

mean_val <- seq(0, 70, length = 1000)
var_th_poi <- var_th_nb <- var_th_qp <- var_th_qp2 <- var_th_qp_exp <- var_th_true <- rep(NA, length(mean_val))
var_th_poi <- mean_val
var_th_nb <- mean_val*(1 + mean_val/phi_nb)
var_th_qp <- psi_qp * mean_val
var_th_qp2 <- psi_qp2 * mean_val^2
var_th_qp_exp <- psi_qp_exp * exp(mean_val)
var_th_true <- psi_true * mean_val

df_obs <- data.frame(mu = mean_obs, v = var_obs, group = as.factor(rep(1, length(mean_obs))))
levels(df_obs$group) <- c("Variance functions")
df_var <- data.frame(mu = mean_val, poi = var_th_poi, nb = var_th_nb,
                     qp = var_th_qp, qp2 = var_th_qp2, qp_exp = var_th_qp_exp)
var_plot <- ggplot(df_obs, aes(x = mu, y = v)) +
  geom_point(shape = 19, size = 2.5) +
  facet_wrap(~ group) +
  geom_line(data = df_var, aes(x = mu, y = poi, color = "poi")) +
  geom_line(data = df_var, aes(x = mu, y = nb, color = "nb")) +
  geom_line(data = df_var, aes(x = mu, y = qp, color = "qp")) +
  geom_line(data = df_var, aes(x = mu, y = qp2, color = "qp2")) +
  geom_line(data = df_var, aes(x = mu, y = qp_exp, color = "qp_exp")) +
  xlim(c(0, 45)) + ylim(c(0, 140)) +
  scale_color_manual(name = "",
                     values = c("poi" = "#339933",
                                "nb" = "#0033CC",
                                "qp" = "#CC3300",
                                "qp2" = "orange",
                                "qp_exp" = "turquoise3"),
                     labels = c(TeX(r'($var(Y) = \mu$)'),
                                TeX(r'($var(Y) = \mu(1 + \mu/\theta)$)'),
                                TeX(r'($var(Y) = \psi\mu$)'),
                                TeX(r'($var(Y) = \psi\mu^2$)'),
                                TeX(r'($var(Y) = \psi e^\mu$)'))) +
  labs(x = "Mean", y = "Variance") +
  theme_light() + 
  theme(legend.text.align = 0,
        strip.text = element_text(size = 16, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"),
        legend.position = "top",
        legend.text = element_text(size = 13))
var_plot

# Error
err <- cbind(c(mean((var_obs - var_poi)^2), mean((var_obs - var_nb)^2),
               mean((var_obs - var_qp)^2), mean((var_obs - var_qp2)^2), mean((var_obs - var_qp_exp)^2)),
             c(mean(abs(var_obs - var_poi)), mean(abs(var_obs - var_nb)),
               mean(abs(var_obs - var_qp)), mean(abs(var_obs - var_qp2)), mean(abs(var_obs - var_qp_exp))))
rownames(err) <- c("poi", "nb", "qp", "qp2", "qp_exp")
colnames(err) <- c("mse", "mae")
err
