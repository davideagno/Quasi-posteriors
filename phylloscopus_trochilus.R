library(lme4)
library(rstan)
library(rstanarm)
library(coda)
library(TeachingDemos)
library(ggplot2)
library(bayesplot)
library(latex2exp)

da <- read.csv("phylloscopus_trochilus.csv")

################################################################################
#                                QUASI-POISSON                                 #
################################################################################

# Preliminary estimates of \psi
hmod <- glmer(y ~ hab + apr_may + year + (1|route), data = da, family = poisson, nAGQ = 0)
beta_mle <- summary(hmod)$coefficients[,1] # EEmpirical risk minimizer
pred_hmod <- predict(hmod, type = "response", re.form = NULL) 
X <- model.matrix(hmod)
stan_data <- as.list(as.data.frame(X))
names(stan_data) <- c("intercept", "hab_Co", "hab_Op", "hab_Urb", "hab_We",
                      "apr_may", "year")
stan_data$y <- da$y
stan_data$N <- N <- NROW(X)
stan_data$p <- P <- NCOL(X)
stan_data$J <- J <- length(unique(da$route))
stan_data$route <- da$route
psi <- sum((da$y - pred_hmod)^2 / pred_hmod) / (stan_data$N - stan_data$p - stan_data$J - 1)
psi # Estimate of dispersion
stan_data$rec_psi <- 1 / psi

burn_in <- 1000
iter <- 5000
fit_ql <- stan(file = 'phylloscopus_trochilus.stan', 
               data = stan_data,
               warmup = burn_in, 
               iter = iter + burn_in, 
               chains = 4, 
               cores = 4)
# Check of convergence
rstan::traceplot(fit_ql, pars = c("beta", "sigma"), inc_warmup = TRUE)
# Save output
HGLM_ql <- as.mcmc(as.matrix(fit_ql))
HGLM_ql_beta <- HGLM_ql[,1:P]
HGLM_ql_delta <- HGLM_ql[,(P+1):(P+J)]
HGLM_ql_sigma <- HGLM_ql[,(P+J+1)]
HGLM_ql_mu <- HGLM_ql[,(P+J+1+N+1):(P+J+1+N+N)]
beta_ql <- colMeans(HGLM_ql_beta)
sigma_ql <- mean(HGLM_ql_sigma)
exp_mu_ql <- colMeans(HGLM_ql_mu)


################################################################################
#                                   POISSON                                    #                        
################################################################################

da_2 <- da; da_2$route <- factor(da_2$route)
fit_poi <- stan_glmer(y ~ hab + apr_may + year + (1|route),
                      data = da_2,
                      family = poisson,
                      warmup = burn_in,
                      iter = iter + burn_in,
                      chains = 4,
                      cores = 4)
# Save outputs
HGLM_poi <- as.mcmc(as.matrix(fit_poi))
HGLM_poi_beta <- HGLM_poi[,1:P]
HGLM_poi_delta <- HGLM_poi[,(P+1):(P+J)]
HGLM_poi_sigma <- HGLM_poi[,(P+J+1)]
beta_poi <- colMeans(HGLM_poi_beta)
sigma_poi <- mean(HGLM_poi_sigma)
exp_mu_poi <- fit_poi$fitted.values


################################################################################
#                              NEGATIVE BINOMIAL                               #                        
################################################################################

fit_nb <- stan_glmer(y ~ hab + apr_may + year + (1|route),
                     data = da_2,
                     family = neg_binomial_2,
                     warmup = burn_in,
                     iter = iter + burn_in,
                     chains = 4,
                     cores = 4)
# Save outputs
HGLM_nb <- as.mcmc(as.matrix(fit_nb))
HGLM_nb_beta <- HGLM_nb[,1:P]
HGLM_nb_delta <- HGLM_nb[,(P+1):(P+J)]
HGLM_nb_sigma <- HGLM_nb[,(P+J+2)]
HGLM_nb_rec_phi <- HGLM_nb[,(P+J+1)]
beta_nb <- colMeans(HGLM_nb_beta) 
sigma_nb <- mean(HGLM_nb_sigma)
rec_phi_nb <- mean(HGLM_nb_rec_phi)
exp_mu_nb <- fit_nb$fitted.values


################################################################################
#                               GOODNESS-OF-FIT                                #                        
################################################################################

mu_nb <- HGLM_nb_beta %*% t(X)
mu_poi <- HGLM_poi_beta %*% t(X)
for (i in 1:NCOL(mu_nb)) {
  mu_nb[,i] <- exp(mu_nb[,i] + HGLM_nb_delta[,da$route[i]])
  mu_poi[,i] <- exp(mu_poi[,i] + HGLM_poi_delta[,da$route[i]])
}
mu_ql_ci <- mu_nb_ci <- mu_poi_ci <- matrix(NA, nrow = N, ncol = 2)
for (i in 1:N) {
  mu_ql_ci[i,] <- emp.hpd(HGLM_ql_mu[,i], conf = 0.95)
  mu_poi_ci[i,] <- emp.hpd(mu_poi[,i], conf = 0.95)
  mu_nb_ci[i,] <- emp.hpd(mu_nb[,i], conf = 0.95)
}
cbind(ql = sum(apply(cbind(mu_ql_ci, da$y), 1, function(x) x[1]<x[3] & x[3]<x[2])) / N,
      poi = sum(apply(cbind(mu_poi_ci, da$y), 1, function(x) x[1]<x[3] & x[3]<x[2])) / N,
      nb = sum(apply(cbind(mu_nb_ci, da$y), 1, function(x) x[1]<x[3] & x[3]<x[2])) / N)


################################################################################
#                                   FIGURE                                     #                        
################################################################################

df_fig <- data.frame(x_ql = c(HGLM_ql_delta[,135],
                              HGLM_ql_mu[,237],
                              HGLM_ql_mu[,236],
                              HGLM_ql_beta[,6],
                              HGLM_ql_sigma),
                     x_poi = c(HGLM_poi_delta[,135],
                               mu_poi[,237],
                               mu_poi[,236],
                               HGLM_poi_beta[,6],
                               HGLM_poi_sigma),
                     x_nb = c(HGLM_nb_delta[,135],
                              mu_nb[,237],
                              mu_nb[,236],
                              HGLM_nb_beta[,6],
                              HGLM_nb_sigma),
                     obs = rep(1:5, each = iter*4),
                     v_line = rep(c(NA, da$y[237], da$y[236], NA, NA), 
                                  each = iter*4))
strip_label = c("Site 135: random intercept", "Site 135, year 2006", 
                "Site 135, year 2008", "Temperature effect", "Random effects variance")
names(strip_label) <- c(1, 2, 3, 4, 5)
fig <- ggplot(df_fig) +
  facet_wrap(~ obs, scales = "free", nrow = 2,
             labeller = labeller(obs = strip_label)) +
  stat_density(aes(x = x_ql, col = "ql"),
               geom = "line", position = "identity", adjust = 2) +
  stat_density(aes(x = x_poi, col = "poi"),
               geom = "line", position = "identity", adjust = 2) +
  stat_density(aes(x = x_nb, col = "nb"),
               geom = "line", position = "identity", adjust = 2) +
  scale_color_manual(name = "Model",
                     values = c("nb" = "#0033CC",
                                "poi" = "#339933",
                                "ql" = "#CC3300"),
                     labels = c("Negative binomial",
                                "Poisson",
                                "Quasi-Poisson")
  ) +
  geom_vline(aes(xintercept = v_line, linetype = "true"), na.rm = TRUE) +
  scale_linetype_manual(name = "",
                        values = c("true" = "dashed"),
                        labels = c("Observed value")) +
  labs(x = "", y = "Density") +
  theme_light() +
  guides(color = guide_legend(order = 1),
         linetype = guide_legend(order = 2)) +
  theme(legend.text.align = 0,
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(size = 12, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"),
        legend.position = "top",
        legend.text = element_text(size = 11))
