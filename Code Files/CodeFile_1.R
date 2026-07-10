# ============================================================================
# COMPLETELY CORRECTED CODE
# Le Cam's Lemma Validation - Two-Sample Rank Tests
# ============================================================================

# Load libraries
library(ggplot2)
library(dplyr)
library(tidyverse)

# ============================================================================
# PART 1: THEORETICAL PARAMETERS
# ============================================================================

# Standardized g0 values
g0_wilcoxon_normal <- (1 / (2 * sqrt(pi))) * sqrt(3)    # 0.4886
g0_wilcoxon_logistic <- (1 / 6) * sqrt(3)               # 0.2887
g0_vdw_normal <- 1                                      # 1.0000
g0_vdw_logistic <- 1 / sqrt(pi)                         # 0.5642

# sigma2_psi
sigma2_psi_normal <- 1
sigma2_psi_logistic <- 1/3

# Function to compute theoretical parameters
compute_params <- function(n, dist, score) {
  
  if (dist == "normal" && score == "wilcoxon") {
    g0 <- g0_wilcoxon_normal
    sigma2_psi <- sigma2_psi_normal
  } else if (dist == "logistic" && score == "wilcoxon") {
    g0 <- g0_wilcoxon_logistic
    sigma2_psi <- sigma2_psi_logistic
  } else if (dist == "normal" && score == "vdw") {
    g0 <- g0_vdw_normal
    sigma2_psi <- sigma2_psi_normal
  } else if (dist == "logistic" && score == "vdw") {
    g0 <- g0_vdw_logistic
    sigma2_psi <- sigma2_psi_logistic
  }
  
  n1 <- n / 2
  n2 <- n / 2
  theta_n <- 1 / sqrt(n)
  
  # Theoretical mean and variance under alternative
  # Note: The 1/sqrt(n) is already in the test statistic!
  mu_n <- theta_n * (n1 * n2 / n) * g0
  sigma2_n <- (n1 * n2 / n) * sigma2_psi
  sigma_n <- sqrt(sigma2_n)
  
  return(list(
    mu_n = mu_n, 
    sigma_n = sigma_n, 
    sigma2_n = sigma2_n,
    theta_n = theta_n, 
    g0 = g0, 
    sigma2_psi = sigma2_psi,
    n1 = n1, 
    n2 = n2,
    n = n,
    dist = dist,
    score = score
  ))
}

# ============================================================================
# PART 2: SIMULATION FUNCTIONS (CORRECTED!)
# ============================================================================

generate_data <- function(n1, n2, dist, theta_n) {
  if (dist == "normal") {
    X <- rnorm(n1, 0, 1)
    Y <- rnorm(n2, theta_n, 1)
  } else if (dist == "logistic") {
    X <- rlogis(n1, 0, 1)
    Y <- rlogis(n2, theta_n, 1)
  }
  return(c(X, Y))
}

compute_T <- function(pooled, n1, n2, score) {
  n <- n1 + n2
  ranks <- rank(pooled, ties.method = "average")
  
  if (score == "wilcoxon") {
    # STANDARDIZED Wilcoxon: sqrt(3) * ranks
    scores <- sqrt(3) * ranks
  } else if (score == "vdw") {
    # Van der Waerden: already standardized
    scores <- qnorm(ranks / (n + 1))
  }
  
  c <- c(rep(0, n1), rep(1, n2))
  c_bar <- mean(c)
  
  # CRITICAL: Divide by sqrt(n) to match asymptotic theory
  T_n <- (1 / sqrt(n)) * sum((c - c_bar) * scores)
  #T_n <- (1 / sqrt(n)) * sum(c  * scores)
  return(T_n)
}

# ============================================================================
# PART 3: RUN SIMULATION
# ============================================================================

run_simulation <- function(n, dist, score, M = 10000, verbose = TRUE) {
  
  params <- compute_params(n, dist, score)
  
  if (verbose) {
    dist_name <- tools::toTitleCase(dist)
    score_name <- ifelse(score == "wilcoxon", "Wilcoxon", "Van der Waerden")
    cat(sprintf("Running: %s, %s, n = %d\n", dist_name, score_name, n))
  }
  
  set.seed(123 + n + 100 * (dist == "logistic") + 1000 * (score == "vdw"))
  
  T_stats <- numeric(M)
  for (r in 1:M) {
    pooled <- generate_data(params$n1, params$n2, dist, params$theta_n)
    T_stats[r] <- compute_T(pooled, params$n1, params$n2, score)
  }
  
  # Empirical moments
  T_mean <- mean(T_stats)
  T_var <- var(T_stats)
  T_sd <- sd(T_stats)
  
  # KS distance
  F_emp <- ecdf(T_stats)
  F_theory <- function(x) pnorm(x, mean = params$mu_n, sd = params$sigma_n)
  ks_values <- sapply(T_stats, function(x) abs(F_emp(x) - F_theory(x)))
  ks <- max(ks_values)
  
  if (verbose) {
    cat(sprintf("  mu_theory = %.4f, mu_emp = %.4f\n", params$mu_n, T_mean))
    cat(sprintf("  sigma_theory = %.4f, sigma_emp = %.4f\n", params$sigma_n, T_sd))
    cat(sprintf("  KS = %.4f\n", ks))
    cat(sprintf("  Var Ratio = %.4f\n", T_var / params$sigma2_n))
  }
  
  return(list(
    T_stats = T_stats,
    mu_n = params$mu_n,
    sigma_n = params$sigma_n,
    sigma2_n = params$sigma2_n,
    T_mean = T_mean,
    T_var = T_var,
    T_sd = T_sd,
    ks = ks,
    n = n,
    dist = dist,
    score = score,
    M = M
  ))
}

# ============================================================================
# PART 4: DENSITY PLOT
# ============================================================================

plot_density <- function(n, dist, score, M = 10000) {
  
  result <- run_simulation(n, dist, score, M, verbose = TRUE)
  
  df <- data.frame(T = result$T_stats)
  
  dist_name <- tools::toTitleCase(result$dist)
  score_name <- ifelse(result$score == "wilcoxon", "Wilcoxon", "Van der Waerden")
  
  p <- ggplot(df, aes(x = T)) +
    geom_density(color = "blue", size = 1.2, kernel = "gaussian") +
    stat_function(
      fun = dnorm, 
      args = list(mean = result$mu_n, sd = result$sigma_n),
      color = "red", size = 1.2, linetype = "dashed"
    ) +
    xlim(-4,4) + 
    ylim(0,1) +
    labs(
      title = paste(dist_name, "Distribution,", score_name, "Score"),
      subtitle = paste0(
        "n = ", result$n,
        "\nKS = ", round(result$ks, 4),
        ", Var Ratio = ", round(result$T_var / result$sigma2_n, 3)
      ),
      x = expression(T[n]),
      y = "Density"
    ) +
    annotate(
      "text", 
      x = -Inf, 
      y = Inf, 
      label = "Blue: Simulated\nRed: Theoretical Normal",
      hjust = -0.1, 
      vjust = 1.5, 
      size = 4
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11),
      axis.title = element_text(size = 12)
    )
  
  return(p)
}

# ============================================================================
# PART 5: COMPARISON PLOT
# ============================================================================
# 
# plot_comparison <- function(n, dist, M = 10000) {
#   
#   dist_name <- tools::toTitleCase(dist)
#   
#   cat("\n========================================\n")
#   cat(dist_name, "COMPARISON at n =", n, "\n")
#   cat("========================================\n")
#   
#   result_w <- run_simulation(n, dist, "wilcoxon", M, verbose = TRUE)
#   result_v <- run_simulation(n, dist, "vdw", M, verbose = TRUE)
#   
#   df_w <- data.frame(T = result_w$T_stats, Score = "Wilcoxon")
#   df_v <- data.frame(T = result_v$T_stats, Score = "Van der Waerden")
#   df_all <- rbind(df_w, df_v)
#   
#   p <- ggplot(df_all, aes(x = T, color = Score, fill = Score)) +
#     geom_density(alpha = 0.3, size = 1.2, kernel = "gaussian") +
#     stat_function(
#       fun = dnorm, 
#       args = list(mean = result_w$mu_n, sd = result_w$sigma_n),
#       color = "blue", size = 0.8, linetype = "dashed"
#     ) +
#     stat_function(
#       fun = dnorm, 
#       args = list(mean = result_v$mu_n, sd = result_v$sigma_n),
#       color = "red", size = 0.8, linetype = "dashed"
#     ) +
#     xlim(-2,2) + 
#     ylim(0,1) +
#     labs(
#       title = paste(dist_name, "Distribution - Comparison at n =", n),
#       subtitle = paste0(
#         "Wilcoxon KS = ", round(result_w$ks, 4), 
#         ", VdW KS = ", round(result_v$ks, 4)
#       ),
#       x = expression(T[n]),
#       y = "Density"
#     ) +
#     scale_color_manual(values = c("Wilcoxon" = "blue", "Van der Waerden" = "red")) +
#     scale_fill_manual(values = c("Wilcoxon" = "blue", "Van der Waerden" = "red")) +
#     theme_minimal() +
#     theme(
#       plot.title = element_text(size = 14, face = "bold"),
#       plot.subtitle = element_text(size = 12),
#       axis.title = element_text(size = 12),
#       legend.position = "bottom"
#     )
#   
#   cat("\n========================================\n")
#   cat("FINAL SUMMARY\n")
#   cat("========================================\n")
#   cat(sprintf("Wilcoxon KS:     %.4f (Var Ratio = %.3f)\n", 
#               result_w$ks, result_w$T_var / result_w$sigma2_n))
#   cat(sprintf("Van der Waerden KS: %.4f (Var Ratio = %.3f)\n", 
#               result_v$ks, result_v$T_var / result_v$sigma2_n))
#   
#   if (result_w$ks < result_v$ks) {
#     cat("\n✓ Wilcoxon performs BETTER for", dist_name, "\n")
#   } else {
#     cat("\n✓ Van der Waerden performs BETTER for", dist_name, "\n")
#   }
#   cat("========================================\n\n")
#   
#   return(p)
# }

# ============================================================================
# RUN THE SIMULATIONS
# ============================================================================

cat("\n")
cat("========================================\n")
cat("RUNNING CORRECTED SIMULATIONS\n")
cat("========================================\n\n")

# Test one scenario first
plot_density(n = 30, dist = "logistic", score = "wilcoxon", M = 5000)
plot_density(n = 30, dist = "logistic", score = "vdw", M = 5000)

plot_density(n = 50, dist = "logistic", score = "wilcoxon", M = 5000)
plot_density(n = 50, dist = "logistic", score = "vdw", M = 5000)




# 
# 
# 
# 
# # Then compare
# plot_comparison(n = 30, dist = "logistic", M = 5000)
# plot_comparison(n = 100, dist = "logistic", M = 5000)
# plot_comparison(n = 30, dist = "normal", M = 5000)
# plot_comparison(n = 100, dist = "normal", M = 5000)
# 
# # ============================================================================
# # END OF SCRIPT
# # ============================================================================