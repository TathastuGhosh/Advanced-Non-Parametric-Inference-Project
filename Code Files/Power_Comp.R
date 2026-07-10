# ============================================================================
# POWER CURVES - VALIDATING THEORY WITH OPTIMAL ASYMPTOTIC CURVES
# Le Cam's Lemma Validation - Two-Sample Rank Tests
# ============================================================================

# Load libraries
library(ggplot2)
library(dplyr)
library(gridExtra)

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

# ============================================================================
# PART 2: THEORETICAL POWER FUNCTION (ASYMPTOTIC ENVELOPE)
# ============================================================================

theoretical_power <- function(theta, g0, sigma2_psi, lambda = 0.5, alpha = 0.05) {
  z_alpha <- qnorm(1 - alpha)
  sigma <- sqrt(lambda * (1 - lambda) * sigma2_psi)
  power <- pnorm(theta * sqrt(lambda * (1 - lambda)) * g0 / sigma - z_alpha)
  return(power)
}

# ============================================================================
# PART 3: SIMULATION FUNCTION
# ============================================================================

get_critical_value <- function(n1, n2, sigma2_psi, alpha = 0.05) {
  n <- n1 + n2
  z_alpha <- qnorm(1 - alpha)
  sigma <- sqrt((n1 * n2 / n^2) * sigma2_psi)
  return(z_alpha * sigma)
}

generate_power_data <- function(n, dist, score, theta_grid, M = 2000, alpha = 0.05) {
  
  # Get theoretical parameters
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
  
  crit_val <- get_critical_value(n1, n2, sigma2_psi, alpha)
  theo_power <- theoretical_power(theta_grid, g0, sigma2_psi, lambda = 0.5, alpha)
  
  emp_power <- numeric(length(theta_grid))
  
  set.seed(123 + n + 100 * (dist == "logistic") + 1000 * (score == "vdw"))
  
  for (j in 1:length(theta_grid)) {
    theta <- theta_grid[j]
    
    reject_count <- 0
    for (r in 1:M) {
      if (dist == "normal") {
        X <- rnorm(n1, 0, 1)
        Y <- rnorm(n2, theta, 1)
      } else {
        X <- rlogis(n1, 0, 1)
        Y <- rlogis(n2, theta, 1)
      }
      pooled <- c(X, Y)
      
      ranks <- rank(pooled, ties.method = "average")
      if (score == "wilcoxon") {
        scores <- sqrt(3) * ranks
      } else {
        scores <- qnorm(ranks / (n + 1))
      }
      
      c <- c(rep(0, n1), rep(1, n2))
      c_bar <- mean(c)
      T_n <- (1 / sqrt(n)) * sum((c - c_bar) * scores)
      
      if (T_n > crit_val) {
        reject_count <- reject_count + 1
      }
    }
    emp_power[j] <- reject_count / M
  }
  
  return(list(
    theta = theta_grid,
    theo_power = theo_power,
    emp_power = emp_power,
    n = n,
    dist = dist,
    score = score
  ))
}

# ============================================================================
# PART 4: RUN SIMULATIONS
# ============================================================================

theta_grid <- seq(0, 3.5, length.out = 25)
sample_sizes <- c(30, 50, 70, 100, 200)
M <- 2000

cat("\n========================================\n")
cat("RUNNING POWER SIMULATIONS\n")
cat("========================================\n")

all_results <- list()
counter <- 1

for (n in sample_sizes) {
  for (dist in c("normal", "logistic")) {
    for (score in c("wilcoxon", "vdw")) {
      cat(sprintf("Running: %s, %s, n = %d\n", dist, score, n))
      result <- generate_power_data(n, dist, score, theta_grid, M = M)
      all_results[[counter]] <- result
      counter <- counter + 1
    }
  }
}

# ============================================================================
# PART 5: CREATE PLOT FOR VAN DER WAERDEN SCORE
# ============================================================================

create_vdw_plot <- function(n_value) {
  
  # Get results for this n
  res_vdw_normal <- all_results[sapply(all_results, function(x) 
    x$n == n_value && x$dist == "normal" && x$score == "vdw")][[1]]
  
  res_vdw_logistic <- all_results[sapply(all_results, function(x) 
    x$n == n_value && x$dist == "logistic" && x$score == "vdw")][[1]]
  
  # Create data frame
  plot_df <- data.frame(
    theta = res_vdw_normal$theta,
    
    # Theoretical asymptotic power for Normal + VdW (Optimal for Normal)
    Theo_Normal_VdW = res_vdw_normal$theo_power,
    
    # Empirical power for Normal + VdW
    Emp_Normal_VdW = res_vdw_normal$emp_power,
    
    # Empirical power for Logistic + VdW (Sub-optimal for Logistic)
    Emp_Logistic_VdW = res_vdw_logistic$emp_power
  )
  
  # Reshape for ggplot
  plot_df_long <- data.frame()
  
  plot_df_long <- rbind(plot_df_long, data.frame(
    theta = plot_df$theta,
    Power = plot_df$Theo_Normal_VdW,
    Curve = "Theoretical: Normal + VdW (Optimal)",
    Type = "Theoretical"
  ))
  
  plot_df_long <- rbind(plot_df_long, data.frame(
    theta = plot_df$theta,
    Power = plot_df$Emp_Normal_VdW,
    Curve = "Empirical: Normal + VdW",
    Type = "Empirical"
  ))
  
  plot_df_long <- rbind(plot_df_long, data.frame(
    theta = plot_df$theta,
    Power = plot_df$Emp_Logistic_VdW,
    Curve = "Empirical: Logistic + VdW",
    Type = "Empirical"
  ))
  
  # Create plot
  p <- ggplot(plot_df_long, aes(x = theta, y = Power, 
                                color = Curve, linetype = Type)) +
    geom_line(size = 1.5) +
    scale_color_manual(
      values = c(
        "Theoretical: Normal + VdW (Optimal)" = "#d62728",
        "Empirical: Normal + VdW" = "#1f77b4",
        "Empirical: Logistic + VdW" = "#2ca02c"
      )
    ) +
    scale_linetype_manual(
      values = c("Theoretical" = "dashed", "Empirical" = "solid"),
      guide = "none"
    ) +
    labs(
      title = paste("Van der Waerden Score - n =", n_value),
      subtitle = "Red dashed: Asymptotic optimal (Normal+VdW) | Blue: Normal+VdW | Green: Logistic+VdW",
      x = expression(theta),
      y = "Power",
      color = "Curve"
    ) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_x_continuous(limits = c(0, 3.5), breaks = seq(0, 3.5, 0.5)) +
    geom_hline(yintercept = 0.05, linetype = "dotted", color = "black", alpha = 0.4) +
    annotate("text", x = 3.2, y = 0.08, label = "alpha = 0.05", size = 3) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 11),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 10),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 9),
      legend.key.width = unit(2, "cm")
    )
  
  return(p)
}

# ============================================================================
# PART 6: CREATE PLOT FOR WILCOXON SCORE
# ============================================================================

create_wilcoxon_plot <- function(n_value) {
  
  # Get results for this n
  res_wilcoxon_logistic <- all_results[sapply(all_results, function(x) 
    x$n == n_value && x$dist == "logistic" && x$score == "wilcoxon")][[1]]
  
  res_wilcoxon_normal <- all_results[sapply(all_results, function(x) 
    x$n == n_value && x$dist == "normal" && x$score == "wilcoxon")][[1]]
  
  # Create data frame
  plot_df <- data.frame(
    theta = res_wilcoxon_logistic$theta,
    
    # Theoretical asymptotic power for Logistic + Wilcoxon (Optimal for Logistic)
    Theo_Logistic_Wilcoxon = res_wilcoxon_logistic$theo_power,
    
    # Empirical power for Logistic + Wilcoxon
    Emp_Logistic_Wilcoxon = res_wilcoxon_logistic$emp_power,
    
    # Empirical power for Normal + Wilcoxon (Sub-optimal for Normal)
    Emp_Normal_Wilcoxon = res_wilcoxon_normal$emp_power
  )
  
  # Reshape for ggplot
  plot_df_long <- data.frame()
  
  plot_df_long <- rbind(plot_df_long, data.frame(
    theta = plot_df$theta,
    Power = plot_df$Theo_Logistic_Wilcoxon,
    Curve = "Theoretical: Logistic + Wilcoxon (Optimal)",
    Type = "Theoretical"
  ))
  
  plot_df_long <- rbind(plot_df_long, data.frame(
    theta = plot_df$theta,
    Power = plot_df$Emp_Logistic_Wilcoxon,
    Curve = "Empirical: Logistic + Wilcoxon",
    Type = "Empirical"
  ))
  
  plot_df_long <- rbind(plot_df_long, data.frame(
    theta = plot_df$theta,
    Power = plot_df$Emp_Normal_Wilcoxon,
    Curve = "Empirical: Normal + Wilcoxon",
    Type = "Empirical"
  ))
  
  # Create plot
  p <- ggplot(plot_df_long, aes(x = theta, y = Power, 
                                color = Curve, linetype = Type)) +
    geom_line(size = 1.5) +
    scale_color_manual(
      values = c(
        "Theoretical: Logistic + Wilcoxon (Optimal)" = "#2ca02c",
        "Empirical: Logistic + Wilcoxon" = "#9467bd",
        "Empirical: Normal + Wilcoxon" = "#1f77b4"
      )
    ) +
    scale_linetype_manual(
      values = c("Theoretical" = "dashed", "Empirical" = "solid"),
      guide = "none"
    ) +
    labs(
      title = paste("Wilcoxon Score - n =", n_value),
      subtitle = "Green dashed: Asymptotic optimal (Logistic+Wilcoxon) | Purple: Logistic+Wilcoxon | Blue: Normal+Wilcoxon",
      x = expression(theta),
      y = "Power",
      color = "Curve"
    ) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_x_continuous(limits = c(0, 3.5), breaks = seq(0, 3.5, 0.5)) +
    geom_hline(yintercept = 0.05, linetype = "dotted", color = "black", alpha = 0.4) +
    annotate("text", x = 3.2, y = 0.08, label = "alpha = 0.05", size = 3) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 11),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 10),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 9),
      legend.key.width = unit(2, "cm")
    )
  
  return(p)
}

# ============================================================================
# PART 7: GENERATE AND DISPLAY ALL PLOTS
# ============================================================================

dir.create("power_plots_validated", showWarnings = FALSE)

for (n in sample_sizes) {
  
  cat(sprintf("\n========================================\n"))
  cat(sprintf("PLOTS FOR n = %d\n", n))
  cat(sprintf("========================================\n"))
  
  # Plot 1: Van der Waerden Score
  p1 <- create_vdw_plot(n)
  print(p1)
  ggsave(sprintf("power_plots_validated/vdw_n%d.png", n), p1, width = 12, height = 8, dpi = 300)
  
  cat(sprintf("Saved: vdw_n%d.png\n", n))
  
  # Plot 2: Wilcoxon Score
  p2 <- create_wilcoxon_plot(n)
  print(p2)
  ggsave(sprintf("power_plots_validated/wilcoxon_n%d.png", n), p2, width = 12, height = 8, dpi = 300)
  
  cat(sprintf("Saved: wilcoxon_n%d.png\n", n))
  
  # Combine both plots side by side
  p_combined <- grid.arrange(
    p1, p2,
    ncol = 2,
    top = grid::textGrob(
      paste("Power Curves - n =", n),
      gp = grid::gpar(fontsize = 18, fontface = "bold")
    )
  )
  
  ggsave(sprintf("power_plots_validated/combined_n%d.png", n), p_combined, width = 16, height = 8, dpi = 300)
  cat(sprintf("Saved: combined_n%d.png\n", n))
  
  cat("\nPress Enter to continue to next n...\n")
  readline()
}

# ============================================================================
# PART 8: SUMMARY
# ============================================================================

cat("\n========================================\n")
cat("ALL PLOTS GENERATED!\n")
cat("========================================\n")
cat("\nFiles saved to 'power_plots_validated/':\n")
cat("\nFor each n = 30, 50, 70, 100, 200:\n")
cat("  - vdw_n[].png        : Van der Waerden score plot\n")
cat("  - wilcoxon_n[].png   : Wilcoxon score plot\n")
cat("  - combined_n[].png   : Both plots side by side\n")
cat("\n========================================\n")

# ============================================================================
# END OF SCRIPT
# ============================================================================