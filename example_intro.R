library(rstan)
library(rstanarm)
library(coda)
library(ggplot2)
library(viridis)

# Simulated dataset
N <- 500
beta_true <- c(1, 1)
psi_true <- 2
set.seed(123)
X <- cbind(rep(1, N),
           rnorm(N, mean = 0, sd = 1))
X <- as.matrix(X)
mu_true = X %*% beta_true
var_true = psi_true * exp(mu_true)
y <- rnorm(N, mean = mu_true, sd = sqrt(var_true))

# Standard Linear Model
fit_lm <- stan_glm(y ~ X - 1,
               family = gaussian,
               prior = normal(0, 10),
               chains = 4,
               cores = 4,
               warmup = 1000,
               iter = 6000)
lm_out <- as.mcmc(as.matrix(fit_lm))
beta_lm <- colMeans(lm_out[,1:2])

# Quasi-posterior
n_logqpost <- function(beta, y, X) {
  mu <- X %*% beta
  - sum((mu - y + 1)*exp(-mu))
}
beta_mle <- nlminb(start = c(0, 0), objective = n_logqpost, y = y, X = X)$par
fit <- X %*% beta_mle
psi <- sum( (y - fit)^2 / exp(fit) ) / (N - 2)
stan_data <- as.list(as.data.frame(X))
names(stan_data) <- c("intercept", "slope")
stan_data$y <- y
stan_data$N <- N
stan_data$p <- 2
stan_data$rec_psi <- 1 / psi
fit_ql <- stan(file = 'example_intro.stan', 
               data = stan_data,
               warmup = 1000, 
               iter = 6000, 
               chains = 4, 
               cores = 4)
rstan::traceplot(fit_ql, pars = c("beta"), inc_warmup = TRUE)
lm_ql <- as.mcmc(as.matrix(fit_ql))
beta_ql <- colMeans(lm_ql[,1:2])

# Figures
df <- data.frame(x = X[,2], y = y, group = rep("Simulated data", length(y)))
df$group <- as.factor(df$group)
plot_1 <- ggplot(data = df, aes(x = x, y = y)) +
  geom_point(size = 1, color = "#336666") + 
  facet_wrap(~ group) +
  labs(x = "X", y = "Y") +
  geom_abline(intercept = beta_lm[1],
              slope = beta_lm[2],
              color = "#0033CC") +
  geom_abline(intercept = beta_ql[1],
              slope = beta_ql[2],
              color = "#CC3300") +
  theme_light() +
  theme(strip.text = element_text(size = 12, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"))

df_ql <- data.frame(cbind(b1 = lm_ql[,1], b2 = lm_ql[,2]),
                    group = rep("Posterior distribution", length(lm_ql[,2])))
df_ql$group <- as.factor(df_ql$group)
df_lm <- data.frame(cbind(b1 = lm_out[,1], b2 = lm_out[,2]))
plot_2 <- ggplot(data = df_ql) +
  facet_wrap(~ group) +
  geom_density_2d(data = df_lm, 
                  aes(x = b1, y = b2, color = "lm"),
                  bins = 4,
                  adjust = 2.5) +
  geom_density_2d(data = df_ql, 
                  aes(x = b1, y = b2, color = "ql"),
                  bins = 4,
                  adjust = 2.5) +
  scale_color_manual(name = "Model",
                     values = c("ql" = "#CC3300",
                                "lm" = "#0033CC"),
                     labels = c("Quasi-posterior", 
                                "Standard posterior")) +
  geom_point(aes(x = beta_true[1], y = beta_true[2], shape = "true"), 
             size = 1.5, stroke = 1.5) + 
  scale_shape_manual(name = "",
                     values = c("true" = 5), 
                     labels = c("True value")) +
  labs(x = "Intercept", y = "Slope") +
  xlim(c(0.65,1.3)) +
  ylim(c(0.55,1.2)) +
  theme_light() +
  theme(strip.text = element_text(size = 12, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93")) +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2))

gridExtra::grid.arrange(plot_1, plot_2, nrow = 1, widths = c(0.42, 0.58))
