################################################################################
#Section 5.3 TABLE 5.3
################################################################################

library(evd)
library(extRemes)

calculate_return_level_formula <- function(T, u, sigma_u, xi, lambda_u, n_y) {
  return_level <- u + (sigma_u / xi) * ((T * n_y * lambda_u)^xi - 1)
  return(return_level)
}

exceedance_probability <- function(threshold, u, sigma_u, xi) {
  if (threshold < u) {
    return(0)
  } else {
    prob <- 1 - evd::pgpd(threshold - u, scale = sigma_u, shape = xi)
    return(prob)
  }
}

annual_exceedance_probability <- function(threshold, u, sigma_u, xi, lambda_u, n_y) {
  prob_exceedance <- exceedance_probability(threshold, u, sigma_u, xi)
  annual_prob <- 1 - (1 - prob_exceedance)^(lambda_u * n_y)
  return(annual_prob * 100)  
}

cumulative_probability <- function(annual_prob, years) {
  cumulative_prob <- 1 - (1 - annual_prob / 100)^years
  return(cumulative_prob * 100)  
}

#We actually do not need CI given we were provided Credible Intervals
#It was for the sake of comparison with other R packages
compute_confidence_intervals <- function(n_samples, u, sigma_u, xi, nreturn, lambda_u, n_y, return_periods) {
  samp <- matrix(rgpd(n_samples * nreturn, u, sigma_u, xi), nreturn, n_samples)
  samp <- apply(samp, 2, sort)
  samp <- apply(samp, 1, sort)
  
  ci_inf <- sapply(return_periods, function(T) {
    return_level_samples <- u + (sigma_u / xi) * ((T * n_y * lambda_u)^xi - 1)
    quantile(return_level_samples, 0.05)
  })
  
  ci_sup <- sapply(return_periods, function(T) {
    return_level_samples <- u + (sigma_u / xi) * ((T * n_y * lambda_u)^xi - 1)
    quantile(return_level_samples, 0.95)
  })
  
  return(list(ci_inf = ci_inf, ci_sup = ci_sup))
}

plot_return_level <- function(data, u, sigma_u, xi, lambda_u, n_y, threshold, main = 'Return Level Plot', xlab = 'Return Period (Years)',
                              ylab = 'Return Level', xlimsup = 100, ci = TRUE, points = TRUE, return_periods = seq(1, 100, length.out = 100), ...) {
  nreturn <- length(data)  
  
  eps <- 10^(-3)
  
  plot(function(T) calculate_return_level_formula(T, u, sigma_u, xi, lambda_u, n_y), from = 1, to = xlimsup, log = 'x', xlab = xlab, ylab = ylab, main = main, ylim = c(100, 1600), ...)
  
  if (points) {
    p_emp <- (1:nreturn - .35) / nreturn
    points(1 / (lambda_u * (1 - p_emp)), sort(data), pch = 1)
  }
  
  if (ci) {
    ci_results <- compute_confidence_intervals(1000, u, sigma_u, xi, nreturn, lambda_u, n_y, return_periods)
    ci_inf <- ci_results$ci_inf
    ci_sup <- ci_results$ci_sup
    lines(return_periods, ci_inf, lty = 2, col = "red")
    lines(return_periods, ci_sup, lty = 2, col = "red")
  }
  
  return_levels_formula <- sapply(return_periods, calculate_return_level_formula, u = u, sigma_u = sigma_u, xi = xi, lambda_u = lambda_u, n_y = n_y)
  lines(return_periods, return_levels_formula, col = "blue", lty = 1)
  
  for (i in seq_along(return_periods)) {
    cat(sprintf("Return Period: %.2f, Return Level: %.2f, Exceedance Probability: %.2f%%\n", 
                return_periods[i], return_levels_formula[i], annual_exceedance_probability(return_levels_formula[i], u, sigma_u, xi, lambda_u, n_y)))
  }
  
  annual_prob <- annual_exceedance_probability(threshold, u, sigma_u, xi, lambda_u, n_y)
  cat(sprintf("Annual probability of exceeding %d nT: %.2f%%\n", threshold, annual_prob))
  
  cumulative_prob_10_years <- cumulative_probability(annual_prob, 10)
  cat(sprintf("Cumulative probability of exceeding %d nT in 10 years: %.2f%%\n", threshold, cumulative_prob_10_years))
}


#Let's recall the estimated parameters from Section 5.1
monitor(as.array(fit_gpd_ds130, pars=c("k","sigma")))

#Ensure this is correctly defined for all 4 datasets 
#4 Datasets : (episode_details$negDst, episode_details$Dual, episode_details_lower$Dual, episode_details_upper$Dual)
#episode_details$negDst, episode_details$Dual, episode_details_lower$Dual will care about 850 nT while episode_details_upper$Dual will care about 1600 nT)

data_rl <- episode_details$negDst  
u <- 130  #Example threshold
sigma_u <- 53.3
xi <- 0.2

#For lambda_u we need [number of declustered exceedances] / [total hourly observations]
#Hence we have 211/582623 = 0.0003621553 (See Coles (2001) or Fawcett for more informations)

lambda_u <- 0.0003621553
n_y <- 365.25 * 24  #Number of hours per year

threshold <- 850  #Threshold to calculate exceedance probability

#Generate return periods from 1 to 100 years
return_periods <- seq(1, 100, length.out = 100)

plot_return_level(data_rl, u, sigma_u, xi, lambda_u, n_y, threshold, return_periods = return_periods)

#Calculate and print the cumulative probabilities of exceeding 850 nT in various years
years <- c(1, 2, 5, 10, 20, 50, 100)
annual_prob <- annual_exceedance_probability(threshold, u, sigma_u, xi, lambda_u, n_y)
cumulative_probs <- sapply(years, cumulative_probability, annual_prob = annual_prob)

cat("Cumulative probabilities of exceeding 850 nT:\n")
for (i in seq_along(years)) {
  cat(sprintf("In %d year(s): %.2f%%\n", years[i], cumulative_probs[i]))
}

#Manually verifying our Return Level
#u + ((53.3 / 0.2) * ((100 * 365.25 * 24 * 0.0003621553)^0.2 - 1))



