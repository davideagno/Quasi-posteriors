library(lme4)
library(rstan)
library(rstanarm)
library(coda)
library(TeachingDemos)
library(ggplot2)
library(bayesplot)
library(latex2exp)
library(plyr)

da <- read.csv("willow_warbler.csv")


################################################################################
#                                QUASI-POISSON                                 #
################################################################################

# Preliminary estimates of \psi
hmod <- glmer(y ~ hab + apr_may + as.factor(year) + (1|route), data = da, family = poisson, nAGQ = 0)
beta_mle <- summary(hmod)$coefficients[,1] # Empirical risk minimizer
pred_hmod <- predict(hmod, type = "response", re.form = NULL) 
X <- model.matrix(hmod)
stan_data <- as.list(as.data.frame(X))
names(stan_data) <- c("intercept", "hab_Co", "hab_Op", "hab_Urb", "hab_We",
                      "apr_may", "year07", "year08")
stan_data$y <- da$y
stan_data$N <- N <- NROW(X)
stan_data$p <- P <- NCOL(X)
stan_data$J <- J <- length(unique(da$route))
stan_data$route <- da$route
psi <- sum((da$y - pred_hmod)^2 / pred_hmod) / (stan_data$N - stan_data$p - stan_data$J - 1)
psi # Estimate of dispersion
stan_data$rec_psi <- 1 / psi
sum(da$y == 0) / N

burn_in <- 500
iter <- 1500
fit_ql <- stan(file = 'willow_warbler.stan', 
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
da_2$year <- factor(da_2$year)
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
phi_nb <- 1/rec_phi_nb
exp_mu_nb <- fit_nb$fitted.values


################################################################################
#                               GOODNESS-OF-FIT                                #                        
################################################################################

mean((da$y - exp_mu_poi)^2 / exp_mu_poi)
mean((da$y - exp_mu_ql)^2 / (psi*exp_mu_ql))
mean((da$y - exp_mu_nb)^2 / (exp_mu_nb*(1 + exp_mu_nb/rec_phi_nb)) )


################################################################################
#                    POSTERIOR MEAN AND CREDIBLE INTERVALS                     #                        
################################################################################

beta_ql_ci <- beta_nb_ci <- beta_poi_ci <- matrix(NA, nrow = P, ncol = 2)
for (i in 1:P) {
  beta_ql_ci[i,] <- emp.hpd(HGLM_ql_beta[,i], conf = 0.95)
  beta_poi_ci[i,] <- emp.hpd(HGLM_poi_beta[,i], conf = 0.95)
  beta_nb_ci[i,] <- emp.hpd(HGLM_nb_beta[,i], conf = 0.95)
}

sigma_ql_ci <- emp.hpd(HGLM_ql_sigma, conf = 0.95)
sigma_poi_ci <- emp.hpd(HGLM_poi_sigma, conf = 0.95)
sigma_nb_ci <- emp.hpd(HGLM_nb_sigma, conf = 0.95)


################################################################################
#                                    FIGURE                                    #                        
################################################################################

idx <- 3
idx_2 <- 59
tot_iter <- NROW(HGLM_ql_beta)
df_delta <- data.frame(x_ql = c(HGLM_ql_delta[,idx], HGLM_ql_delta[,idx_2]),
                       x_poi = c(HGLM_poi_delta[,idx], HGLM_poi_delta[,idx_2]),
                       x_nb = c(HGLM_nb_delta[,idx], HGLM_nb_delta[,idx_2]),
                       obs = c(rep(paste0("(b) Site: ", idx, " - random effect"), tot_iter),
                               rep(paste0("(c) Site: ", idx_2, " - random effect"), tot_iter)))
df_delta$obs <- as.factor(df_delta$obs)
delta_plot <- ggplot(df_delta) +
  facet_wrap(~ obs, scales = "free", nrow = 2) +
  stat_density(aes(x = x_ql, col = "ql"),
               geom = "line", position = "identity", adjust = 2.5) +
  stat_density(aes(x = x_poi, col = "poi"),
               geom = "line", position = "identity", adjust = 2.5) +
  stat_density(aes(x = x_nb, col = "nb"),
               geom = "line", position = "identity", adjust = 2.5) +
  scale_color_manual(name = "Model",
                     values = c("nb" = "#0033CC",
                                "poi" = "#339933",
                                "ql" = "#CC3300"),
                     labels = c("NB",
                                "POI",
                                "QP")
  ) +
  labs(x = TeX(r'($\delta$)'), y = "Density") +
  theme_light() +
  guides(color = guide_legend(order = 1),
         linetype = guide_legend(order = 2)) +
  theme(legend.text.align = 0,
        strip.text = element_text(size = 12, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.title = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top")

sites_loc <- read.csv("sites_location_06_08.csv", head = TRUE)
library(sf)
library(giscoR)
library(dplyr)
fin_map <- gisco_get_countries(resolution = "10", country = "FIN",
                               epsg = "3035") %>% mutate(res = "20M")
df_map <- data.frame(x_site = sites_loc$x_site * 1000,
                     y_site = sites_loc$y_site * 1000,
                     size_site = sites_loc$size_site*2,
                     reff = colMeans(HGLM_ql_delta),
                     hab = sites_loc$hab_site,
                     temp = sites_loc$temp_site,
                     border = rep("gray30", J),
                     lab = rep("(a) Quasi-posterior means for random effects", J))
df_map$lab <- as.factor(df_map$lab)
map_reff <- ggplot(fin_map) +
  geom_sf(fill = "ghostwhite") +
  facet_wrap(~ lab) +
  geom_point(data = df_map, aes(x = x_site, y = y_site, 
                                fill = reff, size = size_site), pch = 23,
             colour = "gray30") +
  geom_point(aes(x = df_map$x_site[idx_2], y = df_map$y_site[idx_2]), pch = 23,
             colour = "red", size = 5, stroke = 1.6) +
  geom_point(aes(x = df_map$x_site[idx], y = df_map$y_site[idx]), pch = 23,
             colour = "blue", size = 5, stroke = 1.6) +
  scale_fill_gradient2(low = "red3", high = "springgreen4", mid = "lightyellow",
                         limits = c(-1, 1)) +
  scale_size_area(breaks = c(2,4,6), label = c(1,2,3)) +
  labs(x = "Longitude", y = "Latitude",
       fill = "Random effect",
       size = "Size") +
  theme_light() +
  theme(strip.text = element_text(size = 12, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 7),
        legend.key.size = unit(11, 'pt'),
        legend.title = element_text(size = 10),
        legend.position = "top")

fig <- gridExtra::grid.arrange(map_reff, delta_plot, ncol = 2)
