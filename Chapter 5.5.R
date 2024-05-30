################################################################################
#Section 5.5
################################################################################

################################################################# GAM diagnostic
library(evgam)
library(evd)
ecdf_function <- ecdf(high_intensity_events$negDst)
ecdf_function(130) #To get a sense for k

fit_gpd_with_best_k <- function(data, cyc_col, excess_col, ks) {
  best_model <- NULL
  best_aic <- Inf
  
  for (k in ks) {
    formula <- list(as.formula(paste(excess_col, "~ s(", cyc_col, ", bs = 'cc', k =", k, ")")), ~ 1)
    tryCatch({
      model <- evgam(formula, data = data, family = "gpd")
      model_aic <- AIC(model)
      
      if (model_aic < best_aic) {
        best_model <- model
        best_aic <- model_aic
      }
    }, error = function(e) {
      message(sprintf("k = %d caused an error: %s", k, e$message))
    })
  }
  
  return(best_model)
}
k_values <- c(1:11)
best_gpd <- fit_gpd_with_best_k(high_intensity_events, "cyc", "excess", k_values)

###################################################################### EVGAM FIT

high_intensity_events <- subset(high_intensity_events, complete.cases(negDst))
zeta <- 0.015 #98 to 99 Quantile, we can choose
theta <- evd::exi(high_intensity_events$Dual, u=70, r = 0)

high_intensity_events_ald <- list(negDst ~ s(cyc, bs = "cc", k = 5), ~ s(cyc, bs = "cc"))
high_ald <- evgam(high_intensity_events_ald, high_intensity_events, family = "ald", ald.args = list(tau = 1 - zeta))
plot(high_ald)
summary(high_ald)
high_intensity_events$threshold <- predict(high_ald)$location
high_intensity_events$excess <- high_intensity_events$negDst - high_intensity_events$threshold
is.na(high_intensity_events$excess[high_intensity_events$excess <= 0]) <- TRUE
summary(high_intensity_events$threshold)

use <- high_intensity_events$year %in% c(1957:2024)
plot(high_intensity_events[use, c("t", "negDst")])
lines(high_intensity_events[use, c("t", "threshold")], col = "red")

high_intensity_events_gpd <- list(excess ~ s(cyc, bs = "cc", k = 5), ~ 1)
high_gpd <- evgam(high_intensity_events_gpd, high_intensity_events, family = "gpd")

rl_df <- data.frame(cyc = seq(0, 365.25)[-1])
rl_df$threshold <- predict(high_ald, rl_df, type = "response")$location
rl_df[, c("psi", "xi")] <- predict(high_gpd, rl_df, type = "response", se.fit = T)
rl_df[, c("psi", "xi")]

FC_sim_ald <- simulate(high_ald, newdata = rl_df, nsim = 1e3, type = "response")
FC_sim_gpd <- simulate(high_gpd, newdata = rl_df, nsim = 1e3, type = "response")
predict(FC_sim_gpd, FC_sim_ald, type = "response", se.fit = F)

rl_100_gpd <- qev(0.99, loc=rl_df$threshold, scale=rl_df$psi, shape=rl_df$xi, m = (365.25*24),  family = "gpd", tau = 1 - zeta)
test <-evir::gpd(high_intensity_events$negDst, threshold = 150)
