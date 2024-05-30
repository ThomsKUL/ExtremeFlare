################################################################################
#Section 4.4
################################################################################

library(ggplot2)
library(gridExtra)
library(poweRlaw)

#Define the Vuong test function without bootstrapping

vuong_test_no_bootstrap <- function(dual_data, negDst_data, threshold = 130) {
  # Fit the models
  lognorm_fit <- conlnorm$new(negDst_data)
  lognorm_fit$setXmin(threshold)
  lognorm_fit$mle()
  
  pareto_fit <- conpl$new(negDst_data)
  pareto_fit$setXmin(threshold)
  pareto_fit$mle()
  
  exp_fit <- conexp$new(negDst_data)
  exp_fit$setXmin(threshold)
  exp_fit$mle()
  
  pareto_truncated_fit <- conpl$new(dual_data)
  pareto_truncated_fit$setXmin(threshold)
  pareto_truncated_fit$mle()
  
  vuong_plt_ln <- compare_distributions(pareto_truncated_fit, lognorm_fit)
  vuong_plt_exp <- compare_distributions(pareto_truncated_fit, exp_fit)
  vuong_pl_ln <- compare_distributions(pareto_fit, lognorm_fit)
  vuong_pl_exp <- compare_distributions(pareto_fit, exp_fit)
  vuong_ln_exp <- compare_distributions(lognorm_fit, exp_fit)
  vuong_plt_pl <- compare_distributions(pareto_truncated_fit, pareto_fit)
  
  results <- list(
    plt_ln_stat = vuong_plt_ln$test_statistic,
    plt_ln_p = vuong_plt_ln$p_two_sided,
    plt_exp_stat = vuong_plt_exp$test_statistic,
    plt_exp_p = vuong_plt_exp$p_two_sided,
    pl_ln_stat = vuong_pl_ln$test_statistic,
    pl_ln_p = vuong_pl_ln$p_two_sided,
    pl_exp_stat = vuong_pl_exp$test_statistic,
    pl_exp_p = vuong_pl_exp$p_two_sided,
    ln_exp_stat = vuong_ln_exp$test_statistic,
    ln_exp_p = vuong_ln_exp$p_two_sided,
    plt_pl_stat = vuong_plt_pl$test_statistic,
    plt_pl_p = vuong_plt_pl$p_two_sided
  )
  
  ratio_plots <- list(
    plt_ln = vuong_plt_ln$ratio,
    plt_exp = vuong_plt_exp$ratio,
    pl_ln = vuong_pl_ln$ratio,
    pl_exp = vuong_pl_exp$ratio,
    ln_exp = vuong_ln_exp$ratio,
    plt_pl = vuong_plt_pl$ratio
  )
  
  return(list(results = results, ratio_plots = ratio_plots))
}


#You just need to switch episode_details_lower with episode_details and episode_details_upper 
#Then plot them adapting their legend

vuong_data <- episode_details_lower
dual_data <- vuong_data$Dual
negDst_data <- vuong_data$negDst
vuong_results <- vuong_test_no_bootstrap(dual_data, negDst_data, threshold = 130)

results <- vuong_results$results
ratio_plots <- vuong_results$ratio_plots

plot_combined_log_likelihood_diff <- function(ratio_plots) {
  combined_data <- data.frame(
    x = c(ratio_plots$plt_ln$x, ratio_plots$plt_exp$x, ratio_plots$pl_ln$x, 
          ratio_plots$pl_exp$x, ratio_plots$ln_exp$x, ratio_plots$plt_pl$x),
    ratio = c(ratio_plots$plt_ln$ratio, ratio_plots$plt_exp$ratio, ratio_plots$pl_ln$ratio, 
              ratio_plots$pl_exp$ratio, ratio_plots$ln_exp$ratio, ratio_plots$plt_pl$ratio),
    comparison = factor(c(
      rep("Truncated Pareto vs. Log-Normal", nrow(ratio_plots$plt_ln)),
      rep("Truncated Pareto vs. Exponential", nrow(ratio_plots$plt_exp)),
      rep("Pareto vs. Log-Normal", nrow(ratio_plots$pl_ln)),
      rep("Pareto vs. Exponential", nrow(ratio_plots$pl_exp)),
      rep("Log-Normal vs. Exponential", nrow(ratio_plots$ln_exp)),
      rep("Truncated Pareto vs. Pareto", nrow(ratio_plots$plt_pl))
    ))
  )
  
  ggplot(combined_data, aes(x = x, y = ratio, color = comparison)) +
    geom_line() +
    ggtitle("Dataset with lower estimates for the 1859 and 2012 events") +
    xlab("|Dst|") +
    ylab("Log-Likelihood Difference") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(color = NULL) +  # This removes the legend title
    scale_color_manual(values = c(
      "Truncated Pareto vs. Log-Normal" = "blue",
      "Truncated Pareto vs. Exponential" = "black",
      "Pareto vs. Log-Normal" = "green",
      "Pareto vs. Exponential" = "red",
      "Log-Normal vs. Exponential" = "purple",
      "Truncated Pareto vs. Pareto" = "orange"
    ))
}

plot_combined_log_likelihood_diff(ratio_plots)
