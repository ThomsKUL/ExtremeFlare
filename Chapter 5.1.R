################################################################################
#Section 5.1
################################################################################

################################################ Functions and other Shenanigans
#See A. Vehtari code for the Stan functions https://github.com/avehtari/BDA_R_demos/blob/master/demos_rstan/gpareto_functions/gpareto.stan
set.seed(100)
library(rstan)
library(loo)
library(dplyr)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(extraDistr)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(rstan)
library(bayesplot)
theme_set(bayesplot::theme_default(base_family = "sans"))

#Ensure the Stan code is saved to a file
stan_code <- '
functions {
  real gpareto_lpdf(vector y, real ymin, real k, real sigma) {
    int N = rows(y);
    real inv_k = inv(k);
    if (k < 0 && max(y - ymin) / sigma > -inv_k) {
      reject("k<0 and max(y-ymin)/sigma > -1/k; found k, sigma =", k, sigma);
    }
    if (sigma <= 0) {
      reject("sigma<=0; found sigma =", sigma);
    }
    if (abs(k) > 1e-15) {
      return -(1 + inv_k) * sum(log1p((y - ymin) * (k / sigma))) - N * log(sigma);
    } else {
      return -sum(y - ymin) / sigma - N * log(sigma);
    }
  }
  real gpareto_cdf(vector y, real ymin, real k, real sigma) {
    real inv_k = inv(k);
    if (k < 0 && max(y - ymin) / sigma > -inv_k) {
      reject("k<0 and max(y-ymin)/sigma > -1/k; found k, sigma =", k, sigma);
    }
    if (sigma <= 0) {
      reject("sigma<=0; found sigma =", sigma);
    }
    if (abs(k) > 1e-15) {
      return exp(sum(log1m_exp((-inv_k) * log1p((y - ymin) * (k / sigma)))));
    } else {
      return exp(sum(log1m_exp(-(y - ymin) / sigma)));
    }
  }
  real gpareto_lcdf(vector y, real ymin, real k, real sigma) {
    real inv_k = inv(k);
    if (k < 0 && max(y - ymin) / sigma > -inv_k) {
      reject("k<0 and max(y-ymin)/sigma > -1/k; found k, sigma =", k, sigma);
    }
    if (sigma <= 0) {
      reject("sigma<=0; found sigma =", sigma);
    }
    if (abs(k) > 1e-15) {
      return sum(log1m_exp((-inv_k) * log1p((y - ymin) * (k / sigma))));
    } else {
      return sum(log1m_exp(-(y - ymin) / sigma));
    }
  }
  real gpareto_lccdf(vector y, real ymin, real k, real sigma) {
    real inv_k = inv(k);
    if (k < 0 && max(y - ymin) / sigma > -inv_k) {
      reject("k<0 and max(y-ymin)/sigma > -1/k; found k, sigma =", k, sigma);
    }
    if (sigma <= 0) {
      reject("sigma<=0; found sigma =", sigma);
    }
    if (abs(k) > 1e-15) {
      return (-inv_k) * sum(log1p((y - ymin) * (k / sigma)));
    } else {
      return -sum(y - ymin) / sigma;
    }
  }
  real gpareto_rng(real ymin, real k, real sigma) {
    if (sigma <= 0) {
      reject("sigma<=0; found sigma =", sigma);
    }
    if (abs(k) > 1e-15) {
      return ymin + (uniform_rng(0, 1) ^ -k - 1) * sigma / k;
    } else {
      return ymin - sigma * log(uniform_rng(0, 1));
    }
  }
}
data {
  real ymin;
  int<lower=0> N;
  vector<lower=ymin>[N] y;
  int<lower=0> Nt;
  vector<lower=ymin>[Nt] yt;
}
transformed data {
  real ymax = max(y);
}
parameters {
  real<lower=0> sigma;
  real<lower=-sigma / (ymax - ymin)> k;
}
model {
  y ~ gpareto(ymin, k, sigma);
}
generated quantities {
  vector[N] log_lik;
  vector[N] yrep;
  vector[Nt] predccdf;
  for (n in 1:N) {
    log_lik[n] = gpareto_lpdf(rep_vector(y[n], 1) | ymin, k, sigma);
    yrep[n] = gpareto_rng(ymin, k, sigma);
  }
  for (nt in 1:Nt) {
    predccdf[nt] = exp(gpareto_lccdf(rep_vector(yt[nt], 1) | ymin, k, sigma));
  }
}
'

#writeLines(stan_code, con = "gpd.stan")
expose_stan_functions("gpd.stan")

###################################################### First Fit + Observed CCDF

#Use this to come back to the original dataset
#high_intensity_events_saved <- high_intensity_events 
#high_intensity_events <- high_intensity_events_saved
high_intensity_events<-episode_details_lower


n <- dim(high_intensity_events)[1]
high_intensity_events <- high_intensity_events %>% arrange(Dual)
high_intensity_events$ccdf <- seq(n,1,-1)/n
quebec_blackout_point <- high_intensity_events[n, ]

ggplot() + 
  geom_point(aes(Dual, ccdf), data = high_intensity_events, size = 1, colour = "blue") +
  coord_trans(x="log10", y="log10", limx=c(100,1000), limy=c(1e-4,1)) +
  scale_y_continuous(breaks=c(1e-5,1e-4,1e-3,1e-2,1e-1,1), limits=c(1e-4,1)) +
  scale_x_continuous(breaks=c(100,200,300,400,500,600,700,850,1000), limits=c(100,1000)) +
  labs(y = "Probability of as Extreme Magnetic Field Disturbance", x= "|Dst|") +
  geom_point(aes(x = quebec_blackout_point$Dual, y = quebec_blackout_point$ccdf), colour = "black", size = 1) +
  geom_text(aes(x = quebec_blackout_point$Dual, y = quebec_blackout_point$ccdf), 
            label = "Quebec blackout", vjust = "top") +
  guides(linetype = FALSE) + 
  theme_bw()

yt130 <- append(10^seq(log10(130), 3, 0.01), 589)
yt130_t <- append(10^seq(log10(130), 3, 0.01), 594)
yt130_tu <- 10^seq(log10(130), log10(1641), length.out = 90)
yt130_tl <- append(10^seq(log10(130), 3, 0.01), 861)

#yt100<-append(10^seq(2,3,.01),850)

ds130<-list(ymin=130, N=dim(episode_details)[1], y=episode_details$negDst, Nt=length(yt130), yt=yt130)
ds130_t<-list(ymin=130, N=dim(episode_details)[1], y=episode_details$Dual, Nt=length(yt130_t), yt=yt130_t)
ds130_tu<-list(ymin=130, N=dim(episode_details_upper)[1], y=episode_details_upper$Dual, Nt=length(yt130_tu), yt=yt130_tu)
ds130_tl<-list(ymin=130, N=dim(episode_details_lower)[1], y=episode_details_lower$Dual, Nt=length(yt130_tl), yt=yt130_tl)

fit_gpd_ds130 <- stan(file='gpd.stan', data=ds130, refresh=0,
                chains=4, seed=100, cores = 4)
str(fit_gpd_ds130)

fit_gpd_ds130_Truncated <- stan(file='gpd.stan', data=ds130_t, refresh=0,
                      chains=4, seed=100, cores = 4)

fit_gpd_ds130_upper <- stan(file='gpd.stan', data=ds130_tu, refresh=0,
                      chains=4, seed=100, cores = 4)

fit_gpd_ds130_lower <- stan(file='gpd.stan', data=ds130_tl, refresh=0,
                      chains=4, seed=100, cores = 4)

monitor(as.array(fit_gpd_ds130, pars=c("k","sigma")))
check_hmc_diagnostics(fit_gpd_ds130)

################################################################################
#Section 5.2
################################################################################

############################################################# Recovery Parameter 
yt130 <- append(10^seq(log10(130), 3, 0.01), 589)
yt130_t <- append(10^seq(log10(130), 3, 0.01), 594)
yt130_tu <- 10^seq(log10(130), log10(1641), length.out = 90)
yt130_tl <- append(10^seq(log10(130), 3, 0.01), 861)


fake_data <- replicate(1e3, gpareto_rng(ymin = 130, 0.2, 53.3))
fake_data_t <- replicate(1e3, gpareto_rng(ymin = 130, 0.2, 53.9))
fake_data_tu <- replicate(1e3, gpareto_rng(ymin = 130, 0.4, 48.8))
fake_data_tl <- replicate(1e3, gpareto_rng(ymin = 130, 0.3, 52.3))

ds_fake<-list(ymin=130, N=1e3, y=fake_data, Nt=length(yt130), yt=yt130)
ds_fake_t<-list(ymin=130, N=1e3, y=fake_data_t, Nt=length(yt130_t), yt=yt130_t)
ds_fake_tu<-list(ymin=130, N=1e3, y=fake_data_tu, Nt=length(yt130_tu), yt=yt130_tu)
ds_fake_tl<-list(ymin=130, N=1e3, y=fake_data_tl, Nt=length(yt130_tl), yt=yt130_tl)

fake_gpd_ds130 <- stan(file='gpd.stan', data=ds_fake, refresh=0,
                      chains=4, seed=100, cores = 4)

fake_gpd_ds130_Truncated <- stan(file='gpd.stan', data=ds_fake_t, refresh=0,
                                chains=4, seed=100, cores = 4)

fake_gpd_ds130_upper <- stan(file='gpd.stan', data=ds_fake_tu, refresh=0,
                            chains=4, seed=100, cores = 4)

fake_gpd_ds130_lower <- stan(file='gpd.stan', data=ds_fake_tl, refresh=0,
                            chains=4, seed=100, cores = 4)

monitor(as.array(fake_gpd_ds130_lower, pars=c("k","sigma")))
check_hmc_diagnostics(fake_gpd_ds130)

######################################## Plot the Recovery (Figures 5.1 and 5.2)

posterior_classical <- as.matrix(fake_gpd_ds130, pars = c("k", "sigma"))
posterior_truncated <- as.matrix(fake_gpd_ds130_Truncated, pars = c("k", "sigma"))
posterior_upper <- as.matrix(fake_gpd_ds130_upper, pars = c("k", "sigma"))
posterior_lower <- as.matrix(fake_gpd_ds130_lower, pars = c("k", "sigma"))

posterior_classical <- data.frame(posterior_classical, Distribution = "Classical GPD")
posterior_truncated <- data.frame(posterior_truncated, Distribution = "Truncated GPD")
posterior_upper <- data.frame(posterior_upper, Distribution = "Truncated GPD+Higher estimates")
posterior_lower <- data.frame(posterior_lower, Distribution = "Truncated GPD+Lower estimates")

posterior_all<-NULL
posterior_all <- rbind(posterior_classical, posterior_truncated, posterior_upper, posterior_lower)
colnames(posterior_all) <- c("k", "sigma", "Distribution")

true_values <- data.frame(
  Distribution = c("Classical GPD", "Truncated GPD", "Truncated GPD+Higher estimates", "Truncated GPD+Lower estimates"),
  k = c(0.20, 0.20, 0.40, 0.30),
  sigma = c(53.30, 53.90, 48.20, 52.00)
)

plot_k <- ggplot(posterior_all, aes(x = k)) +
  geom_histogram(fill = "blue", alpha = 0.6, position = "identity", bins = 30) +
  geom_vline(data = true_values, aes(xintercept = k), color = "red", linetype = "dashed") +
  labs(title = "Recovery of the Shape parameter",
       x = "k", y = "Frequency") +
  theme_minimal() +
  facet_wrap(~ Distribution, scales = "free_y") +
  theme(legend.position = "none")

plot_sigma <- ggplot(posterior_all, aes(x = sigma)) +
  geom_histogram(fill = "blue", alpha = 0.6, position = "identity", bins = 30) +
  geom_vline(data = true_values, aes(xintercept = sigma), color = "red", linetype = "dashed") +
  labs(title = "Recovery of the Scale parameter",
       x = "sigma", y = "Frequency") +
  theme_minimal() +
  facet_wrap(~ Distribution, scales = "free_y") +
  theme(legend.position = "none")

print(plot_k)
print(plot_sigma)


############# Plot Model ability to predict extreme values (Figures 5.3 and 5.4)

create_ppc2_plot <- function(fit, observed_data, title) {
  gpd_params <- rstan::extract(fit)
  max_log_dst <- max(log(observed_data$Dual))
  max_log_replicates <- apply(log(gpd_params$yrep), 1, max)
  
  max_log_data <- data.frame(
    Observed = max_log_dst,
    Replicates = max_log_replicates
  )
  
  ppc2_plot <- ggplot(max_log_data, aes(x = Replicates)) +
    geom_histogram(fill = "blue", alpha = 0.6, bins = 30) +
    geom_vline(aes(xintercept = Observed), color = "red", linetype = "dashed") +
    labs(x = "Max Log Magnitude", y = "Frequency", title = title) +
    theme_minimal()
  
  return(ppc2_plot)
}

ppc2_classical_lower <- create_ppc2_plot(fit_gpd_ds130, episode_details_lower, "Classical GPD")
ppc2_truncated_lower <- create_ppc2_plot(fit_gpd_ds130_Truncated, episode_details_lower, "Truncated GPD")
ppc2_upper_lower <- create_ppc2_plot(fit_gpd_ds130_upper, episode_details_lower, "Truncated GPD+Higher estimates")
ppc2_lower_lower <- create_ppc2_plot(fit_gpd_ds130_lower, episode_details_lower, "Truncated GPD+Lower estimates")

grid_arrange_lower <- grid.arrange(ppc2_classical_lower, ppc2_truncated_lower, ppc2_upper_lower, ppc2_lower_lower, ncol = 2)

ppc2_classical_upper <- create_ppc2_plot(fit_gpd_ds130, episode_details_upper, "Classical GPD")
ppc2_truncated_upper <- create_ppc2_plot(fit_gpd_ds130_Truncated, episode_details_upper, "Truncated GPD")
ppc2_upper_upper <- create_ppc2_plot(fit_gpd_ds130_upper, episode_details_upper, "Truncated GPD+Higher estimates")
ppc2_lower_upper <- create_ppc2_plot(fit_gpd_ds130_lower, episode_details_upper, "Truncated GPD+Lower estimates")

grid_arrange_upper <- grid.arrange(ppc2_classical_upper, ppc2_truncated_upper, ppc2_upper_upper, ppc2_lower_upper, ncol = 2)

########################################################### Goodness of fit test
#Reitere with each datasets 
#Datasets: (episode_details$Dual, episode_details_lower$Dual, episode_details_upper$Dual)
library(eva)
eva::gpdAd(
  episode_details$Dual,
  bootstrap = TRUE,
  bootnum = 10000,
  allowParallel = FALSE,
  numCores = 1
)

eva::gpdCvm(
  episode_details$Dual,
  bootstrap = TRUE,
  bootnum = 10000,
  allowParallel = FALSE,
  numCores = 1
)


################################################################################
#Section 5.3
################################################################################

######################################################## Plot the estimated CCDF
high_intensity_events<-episode_details_upper

n <- dim(high_intensity_events)[1]
high_intensity_events <- high_intensity_events %>% arrange(Dual)
high_intensity_events$ccdf <- seq(n,1,-1)/n
quebec_blackout_point <- high_intensity_events[n, ]

extract_predccdf <- function(fit, yt, label) {
  gpd_params <- rstan::extract(fit)
  num_rows <- ncol(gpd_params$predccdf)
  yt_trimmed <- yt[1:num_rows]
  mu <- apply(t(gpd_params$predccdf), 1, quantile, c(0.05, 0.5, 0.95)) %>%
    t() %>% data.frame(x = yt_trimmed, .)
  mu$fit <- label
  return(mu)
}

# Extract predicted CCDFs for each fit
yt130 <- append(10^seq(log10(130), 3, 0.01), 589)
yt130_t <- append(10^seq(log10(130), 3, 0.01), 594)
yt130_tu <- 10^seq(log10(130), log10(1641), length.out = 90)
yt130_tl <- append(10^seq(log10(130), 3, 0.01), 861)


mu_gpd <- extract_predccdf(fake_gpd_ds130, yt130, "Classical GPD")
mu_truncated <- extract_predccdf(fake_gpd_ds130_Truncated, yt130_t, "Truncated GPD")
mu_upper <- extract_predccdf(fake_gpd_ds130_upper, yt130_tu, "Truncated GPD with added upper obs")
mu_lower <- extract_predccdf(fake_gpd_ds130_lower, yt130_tl, "Truncated GPD with added lower obs")

mu_gpd <- extract_predccdf(fit_gpd_ds130, yt130, "Classical GPD")
mu_truncated <- extract_predccdf(fit_gpd_ds130_Truncated, yt130_t, "Truncated GPD")
mu_upper <- extract_predccdf(fit_gpd_ds130_upper, yt130_tu, "Truncated GPD with added upper obs")
mu_lower <- extract_predccdf(fit_gpd_ds130_lower, yt130_tl, "Truncated GPD with added lower obs")

mu_combined<-NULL
mu_combined <- bind_rows(mu_gpd, mu_truncated, mu_upper, mu_lower) %>%
  gather(pct, y, -x, -fit)

fit_colors <- c("Classical GPD" = "green", "Truncated GPD" = "orange", "Truncated GPD with added upper obs" = "purple", "Truncated GPD with added lower obs" = "red")

ggplot() + 
  geom_point(aes(Dual, ccdf), data = high_intensity_events, size = 1.2, color = "black") +  
  geom_line(aes(x, y, linetype = pct, color = fit), data = mu_combined, size = 0.7) +  # Increased line thickness
  scale_linetype_manual(values = c(2, 1, 2)) +
  scale_color_manual(values = fit_colors) +
  coord_trans(x = "log10", y = "log10", limx = c(100, 1800), limy = c(1e-4, 1)) +  # Adjust x-axis limits to 100-1800
  scale_y_continuous(breaks = c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1), limits = c(1e-4, 1)) +
  scale_x_continuous(breaks = c(100, 200, 300, 400, 500, 600, 700, 850, 1000, 1182, 1400, 1600, 1800), limits = c(100, 1800)) +
  geom_text(aes(x = quebec_blackout_point$Dual, y = quebec_blackout_point$ccdf * 0.95), label = "1859 Carrington", vjust = "top", nudge_y = -0.0005) +
  geom_text(aes(x = high_intensity_events[n-1, ]$Dual, y = high_intensity_events[n-1, ]$ccdf * 0.95), label = "2012 Near Miss", vjust = "top", nudge_y = -0.0005) +
  geom_text(aes(x = high_intensity_events[n-2, ]$Dual, y = high_intensity_events[n-2, ]$ccdf * 0.95), label = "1989 Quebec", vjust = "top", nudge_y = -0.0005) +
  labs(y = "Probability of as Extreme Magnetic Field Disturbance", x = "|Dst|", color = "Fit") +
  guides(linetype = FALSE) + 
  theme_bw() +
  theme(legend.position = "bottom")

################################################################################
#Section 5.4
################################################################################

################################################################ Noisy estimates
high_intensity_events<-episode_details


yt130 <- append(10^seq(log10(130), 3, 0.01), 589)
yt130_t <- append(10^seq(log10(130), 3, 0.01), 594)
yt130_tu <- 10^seq(log10(130), log10(1641), length.out = 90)
yt130_tl <- append(10^seq(log10(130), 3, 0.01), 861)

dsn130<-list(ymin=130, N=dim(episode_details)[1], y=episode_details$negDst, Nt=length(yt130), yt=yt130, noise=0.2*high_intensity_events$Dual)
dsn130_t<-list(ymin=130, N=dim(episode_details)[1], y=episode_details$Dual, Nt=length(yt130_t), yt=yt130_t, noise=0.2*high_intensity_events$Dual)
dsn130_tu<-list(ymin=130, N=dim(episode_details_upper)[1], y=episode_details_upper$Dual, Nt=length(yt130_tu), yt=yt130_tu, noise=0.2*high_intensity_events$Dual)
dsn130_tl<-list(ymin=130, N=dim(episode_details_lower)[1], y=episode_details_lower$Dual, Nt=length(yt130_tl), yt=yt130_tl, noise=0.2*high_intensity_events$Dual)


fit_gpd_ds130_noise <- stan(file='gpd.stan', data=dsn130, refresh=0,
                      chains=4, seed=100, cores = 4, iter = 5000)

fit_gpd_ds130_Truncated_noise <- stan(file='gpd.stan', data=dsn130_t, refresh=0,
                                chains=4, seed=100, cores = 4, iter = 5000)

fit_gpd_ds130_upper_noise <- stan(file='gpd.stan', data=dsn130_tu, refresh=0,
                            chains=4, seed=100, cores = 4, iter = 5000)

fit_gpd_ds130_lower_noise <- stan(file='gpd.stan', data=dsn130_tl, refresh=0,
                            chains=4, seed=100, cores = 4, iter = 5000)

monitor(as.array(fit_gpd_ds130_noise, pars=c("k","sigma")))
monitor(as.array(fit_gpd_ds130_noise, pars=c("k","sigma")))
monitor(as.array(fit_gpd_ds130_noise, pars=c("k","sigma")))


##################################################################### Plot noise
library(rstan)
library(bayesplot)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(grid)

extract_posterior_samples <- function(fit, fit_name) {
  posterior_samples <- as.matrix(fit, pars = c("k", "sigma"))
  posterior_df <- as.data.frame(posterior_samples)
  posterior_df$fit <- fit_name
  return(posterior_df)
}

posterior_classical <- extract_posterior_samples(fit_gpd_ds130, "Classical")
posterior_classical_noisy <- extract_posterior_samples(fit_gpd_ds130_noise, "Noisy Classical")
compare_classical <- bind_rows(posterior_classical, posterior_classical_noisy)

posterior_truncated <- extract_posterior_samples(fit_gpd_ds130_Truncated, "Truncated")
posterior_truncated_noisy <- extract_posterior_samples(fit_gpd_ds130_Truncated_noise, "Noisy Truncated")
compare_truncated <- bind_rows(posterior_truncated, posterior_truncated_noisy)

posterior_upper <- extract_posterior_samples(fit_gpd_ds130_upper, "Upper")
posterior_upper_noisy <- extract_posterior_samples(fit_gpd_ds130_upper_noise, "Noisy Upper")
compare_upper <- bind_rows(posterior_upper, posterior_upper_noisy)

posterior_lower <- extract_posterior_samples(fit_gpd_ds130_lower, "Lower")
posterior_lower_noisy <- extract_posterior_samples(fit_gpd_ds130_lower_noise, "Noisy Lower")
compare_lower <- bind_rows(posterior_lower, posterior_lower_noisy)

prepare_data_for_plotting <- function(compare_data) {
  compare_data <- compare_data %>%
    pivot_longer(cols = c(k, sigma), names_to = "parameter") %>%
    mutate(parameter = recode(parameter, k = "xi"))
  return(compare_data)
}

compare_classical_long <- prepare_data_for_plotting(compare_classical)
compare_truncated_long <- prepare_data_for_plotting(compare_truncated)
compare_upper_long <- prepare_data_for_plotting(compare_upper)
compare_lower_long <- prepare_data_for_plotting(compare_lower)

legend_labels <- c("Classical" = "Classical GPD",
                   "Noisy Classical" = "Noisy GPD Classical",
                   "Truncated" = "Truncated GPD",
                   "Noisy Truncated" = "Noisy Truncated GPD",
                   "Upper" = "Truncated GPD+Upper estimates",
                   "Noisy Upper" = "Truncated Noisy GPD+Upper estimates",
                   "Lower" = "Truncated GPD+Lower estimates",
                   "Noisy Lower" = "Noisy Truncated GPD+Lower estimates")

plot_density <- function(data, legend_labels) {
  ggplot(data, aes(x = value, color = fit, fill = fit)) +
    geom_density(alpha = 0.3) +  # Adjust transparency with alpha
    facet_wrap(~parameter, scales = "free", labeller = labeller(parameter = label_parsed)) +  # Use labeller to parse expressions
    labs(x = "Parameter Value", y = expression(f(x))) +  # Use mathematical notation for y-axis label
    scale_color_manual(values = c("Classical" = "blue", "Noisy Classical" = "red",
                                  "Truncated" = "green", "Noisy Truncated" = "purple",
                                  "Upper" = "orange", "Noisy Upper" = "brown",
                                  "Lower" = "pink", "Noisy Lower" = "grey"),
                       labels = legend_labels) +
    scale_fill_manual(values = c("Classical" = "blue", "Noisy Classical" = "red",
                                 "Truncated" = "green", "Noisy Truncated" = "purple",
                                 "Upper" = "orange", "Noisy Upper" = "brown",
                                 "Lower" = "pink", "Noisy Lower" = "grey"),
                      labels = legend_labels) +
    theme_minimal() +  # Use a clean and academic theme
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),  # Remove legend title
      legend.text = element_text(size = 10),  # Reduce legend text size
      legend.box.spacing = unit(0.5, "lines"),  # Reduce spacing between legend items
      strip.text = element_text(face = "bold", size = 14),  # Bold and increase facet titles size
      panel.grid.major = element_line(size = 0.5, linetype = 'solid', color = 'grey80'),  # Major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      panel.spacing = unit(1, "lines")  # Space between panels
    ) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE), fill = guide_legend(nrow = 2, byrow = TRUE))  # Arrange legend in two rows
}

plot1 <- plot_density(compare_classical_long, legend_labels)
plot2 <- plot_density(compare_truncated_long, legend_labels)
plot3 <- plot_density(compare_upper_long, legend_labels)
plot4 <- plot_density(compare_lower_long, legend_labels)

grid.arrange(
  arrangeGrob(plot1 + labs(title = NULL), plot2 + labs(title = NULL), plot3 + labs(title = NULL), plot4 + labs(title = NULL), ncol = 2),
  top = textGrob("Posterior Densities with Gaussian noise", gp = gpar(fontface = "bold", fontsize = 16))
)

################################################ Plot Pareto Smoothed Importance
min(1-(1/log10(211)),0.7)

gpd_params1 <- rstan::extract(fit_gpd_ds130)
gpd_params2 <- rstan::extract(fit_gpd_ds130_Truncated)
gpd_params3 <- rstan::extract(fit_gpd_ds130_upper)
gpd_params4 <- rstan::extract(fit_gpd_ds130_lower)

psiss1 <- psis(-gpd_params1$log_lik)
psiss2 <- psis(-gpd_params2$log_lik)
psiss3 <- psis(-gpd_params3$log_lik)
psiss4 <- psis(-gpd_params4$log_lik)

clrs <- color_scheme_get("brightblue")

n1 <- length(psiss1$diagnostics$pareto_k)
n2 <- length(psiss2$diagnostics$pareto_k)
n3 <- length(psiss3$diagnostics$pareto_k)
n4 <- length(psiss4$diagnostics$pareto_k)

pkhats <- ggplot() + 
  geom_point(aes(x = seq(1, n1), y = psiss1$diagnostics$pareto_k), color = clrs[[5]]) + 
  labs(y = "Pareto"~hat(xi), x = "Classical GPD") + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  ylim(-0.25, 1) + 
  theme_default()

pkhats <- ggplot() + 
  geom_point(aes(x = seq(1, n1), y = psiss1$diagnostics$pareto_k), color = clrs[[5]]) + 
  labs(y = expression("Pareto"~hat(xi)), x = "Classical GPD") + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_hline(yintercept = 0.57, linetype = "dotted", color = "red") +
  annotate("text", x = n1 * 0.98, y = 0.62, label = expression("Max"~hat(k)~"for reliable PSIS"), color = "red", hjust = 1) +
  ylim(-0.25, 1) + 
  theme_default()

pkhats2 <- ggplot() + 
  geom_point(aes(x = seq(1, n2), y = psiss2$diagnostics$pareto_k), color = clrs[[5]]) + 
  labs(y = expression("Pareto"~hat(xi)), x = "Truncated GPD") + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_hline(yintercept = 0.57, linetype = "dotted", color = "red") +
  ylim(-0.25, 1) + 
  theme_default()

pkhats3 <- ggplot() + 
  geom_point(aes(x = seq(1, n3), y = psiss3$diagnostics$pareto_k), color = clrs[[5]]) + 
  labs(y = expression("Pareto"~hat(xi)), x = "Truncated GPD+Upper estimates") + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_hline(yintercept = 0.57, linetype = "dotted", color = "red") +
  ylim(-0.25, 1) + 
  theme_default()

pkhats4 <- ggplot() + 
  geom_point(aes(x = seq(1, n4), y = psiss4$diagnostics$pareto_k), color = clrs[[5]]) + 
  labs(y = expression("Pareto"~hat(xi)), x = "Truncated GPD+Lower estimates") + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_hline(yintercept = 0.57, linetype = "dotted", color = "red") +
  ylim(-0.25, 1) + 
  theme_default()

grid.arrange(pkhats, pkhats2, pkhats3, pkhats4, ncol = 2)

################################################################################
#Appendix
################################################################################

#################################################################### Convergence

color_scheme_set("viridis")

posterior1 <- as.array(fit_gpd_ds130)
posterior2 <- as.array(fit_gpd_ds130_Truncated)
posterior3 <- as.array(fit_gpd_ds130_upper)
posterior4 <- as.array(fit_gpd_ds130_lower)

chain1 <- mcmc_rank_ecdf(posterior1, prob = 0.99, pars = c("k", "sigma")) + labs(title = "Classical GPD")
chain2 <- mcmc_rank_ecdf(posterior2, prob = 0.99, pars = c("k", "sigma")) + labs(title = "Truncated GPD")
chain3 <- mcmc_rank_ecdf(posterior3, prob = 0.99, pars = c("k", "sigma")) + labs(title = "Truncated GPD with added upper obs")
chain4 <- mcmc_rank_ecdf(posterior4, prob = 0.99, pars = c("k", "sigma")) + labs(title = "Truncated GPD with added lower obs")

grid_arrange_upper <- grid.arrange(chain1 + theme(axis.title.x = element_blank(), axis.title.y = element_blank()),
                                   chain2 + theme(axis.title.x = element_blank(), axis.title.y = element_blank()),
                                   chain3 + theme(axis.title.x = element_blank(), axis.title.y = element_blank()),
                                   chain4 + theme(axis.title.x = element_blank(), axis.title.y = element_blank()),
                                   ncol = 2)

print(grid_arrange_upper)


########################################################################## Shiny
#install.packages("shinystan")
library(shinystan)
launch_shinystan(fit_gpd_ds130)

################################################################################
