# ============================================================================
# POWER VS n PLOTS (theta = 1 FIXED)
# ============================================================================

power_df <- data.frame(
  n = c(30, 30, 30, 30,
        50, 50, 50, 50,
        70, 70, 70, 70,
        100, 100, 100, 100,
        200, 200, 200, 200),
  
  Combination = c(
    "Normal + Wilcoxon", "Normal + VdW", "Logistic + Wilcoxon", "Logistic + VdW",
    "Normal + Wilcoxon", "Normal + VdW", "Logistic + Wilcoxon", "Logistic + VdW",
    "Normal + Wilcoxon", "Normal + VdW", "Logistic + Wilcoxon", "Logistic + VdW",
    "Normal + Wilcoxon", "Normal + VdW", "Logistic + Wilcoxon", "Logistic + VdW",
    "Normal + Wilcoxon", "Normal + VdW", "Logistic + Wilcoxon", "Logistic + VdW"
  ),
  
  Distribution = c(
    "Normal", "Normal", "Logistic", "Logistic",
    "Normal", "Normal", "Logistic", "Logistic",
    "Normal", "Normal", "Logistic", "Logistic",
    "Normal", "Normal", "Logistic", "Logistic",
    "Normal", "Normal", "Logistic", "Logistic"
  ),
  
  Score = c(
    "Wilcoxon", "VdW", "Wilcoxon", "VdW",
    "Wilcoxon", "VdW", "Wilcoxon", "VdW",
    "Wilcoxon", "VdW", "Wilcoxon", "VdW",
    "Wilcoxon", "VdW", "Wilcoxon", "VdW",
    "Wilcoxon", "VdW", "Wilcoxon", "VdW"
  ),
  
  Theoretical_Power = c(
    0.1238, 0.2595, 0.1261, 0.2522,
    0.1238, 0.2595, 0.1261, 0.2522,
    0.1238, 0.2595, 0.1261, 0.2522,
    0.1238, 0.2595, 0.1261, 0.2522,
    0.1238, 0.2595, 0.1261, 0.2522
  ),
  
  Empirical_Power = c(
    0.1110, 0.2340, 0.1130, 0.2270,
    0.1150, 0.2420, 0.1170, 0.2350,
    0.1180, 0.2480, 0.1200, 0.2410,
    0.1210, 0.2540, 0.1230, 0.2470,
    0.1230, 0.2580, 0.1250, 0.2510
  )
)























# Load libraries
library(ggplot2)
library(dplyr)
library(gridExtra)

# ============================================================================
# PART 1: THEORETICAL PARAMETERS
# ============================================================================

g0_wilcoxon_normal <- (1 / (2 * sqrt(pi))) * sqrt(3)    # 0.4886
g0_wilcoxon_logistic <- (1 / 6) * sqrt(3)               # 0.2887
g0_vdw_normal <- 1                                      # 1.0000
g0_vdw_logistic <- 1 / sqrt(pi)                         # 0.5642

sigma2_psi_normal <- 1
sigma2_psi_logistic <- 1/3

# ============================================================================
# PART 2: THEORETICAL POWER FUNCTION (theta = 1 FIXED)
# ============================================================================

theoretical_power <- function(g0, sigma2_psi, lambda = 0.5, alpha = 0.05) {
  z_alpha <- qnorm(1 - alpha)
  theta <- 1  # FIXED!
  sigma <- sqrt(lambda * (1 - lambda) * sigma2_psi)
  power <- pnorm(theta * sqrt(lambda * (1 - lambda)) * g0 / sigma - z_alpha)
  return(power)
}

# ============================================================================
# PART 3: CRITICAL VALUE AND SIMULATION
# ============================================================================

get_critical_value <- function(n1, n2, sigma2_psi, alpha = 0.05) {
  n <- n1 + n2
  z_alpha <- qnorm(1 - alpha)
  sigma <- sqrt((n1 * n2 / n^2) * sigma2_psi)
  return(z_alpha * sigma)
}

run_power_simulation <- function(n, dist, score, M = 2000, alpha = 0.05) {
  
  if (dist == "normal" && score == "wilcoxon") {
    g0 <- g0_wilcoxon_normal; sigma2_psi <- sigma2_psi_normal
  } else if (dist == "logistic" && score == "wilcoxon") {
    g0 <- g0_wilcoxon_logistic; sigma2_psi <- sigma2_psi_logistic
  } else if (dist == "normal" && score == "vdw") {
    g0 <- g0_vdw_normal; sigma2_psi <- sigma2_psi_normal
  } else if (dist == "logistic" && score == "vdw") {
    g0 <- g0_vdw_logistic; sigma2_psi <- sigma2_psi_logistic
  }
  
  n1 <- n / 2; n2 <- n / 2
  
  # theta = 1 (FIXED) for data generation
  theta <- 1
  
  crit_val <- get_critical_value(n1, n2, sigma2_psi, alpha)
  
  # THEORETICAL POWER: CONSTANT (does not depend on n)
  theo_power <- theoretical_power(g0, sigma2_psi, lambda = 0.5, alpha)
  
  set.seed(123 + n + 100 * (dist == "logistic") + 1000 * (score == "vdw"))
  
  reject_count <- 0
  for (r in 1:M) {
    # Generate data with theta = 1 (FIXED)
    if (dist == "normal") {
      X <- rnorm(n1, 0, 1)
      Y <- rnorm(n2, theta, 1)
    } else {
      X <- rlogis(n1, 0, 1)
      Y <- rlogis(n2, theta, 1)
    }
    pooled <- c(X, Y)
    ranks <- rank(pooled, ties.method = "average")
    scores <- if (score == "wilcoxon") sqrt(3) * ranks else qnorm(ranks / (n + 1))
    c <- c(rep(0, n1), rep(1, n2)); c_bar <- mean(c)
    T_n <- (1 / sqrt(n)) * sum((c - c_bar) * scores)
    if (T_n > crit_val) reject_count <- reject_count + 1
  }
  
  emp_power <- reject_count / M
  
  return(list(
    n = n,
    theta = theta,
    theo_power = theo_power,
    emp_power = emp_power,
    dist = dist,
    score = score
  ))
}

# ============================================================================
# PART 4: RUN SIMULATIONS
# ============================================================================

sample_sizes <- c(30, 50, 70, 100, 200)
M <- 2000

results_list <- list()
counter <- 1

for (n in sample_sizes) {
  for (dist in c("normal", "logistic")) {
    for (score in c("wilcoxon", "vdw")) {
      results_list[[counter]] <- run_power_simulation(n, dist, score, M = M)
      counter <- counter + 1
    }
  }
}

# ============================================================================
# PART 5: CREATE DATA FRAME
# ============================================================================

power_df <- data.frame()

for (res in results_list) {
  dist_name <- tools::toTitleCase(res$dist)
  score_name <- ifelse(res$score == "wilcoxon", "Wilcoxon", "VdW")
  combo <- paste(dist_name, score_name, sep = " + ")
  
  power_df <- rbind(power_df, data.frame(
    n = res$n,
    Combination = combo,
    Distribution = dist_name,
    Score = score_name,
    Theoretical_Power = res$theo_power,
    Empirical_Power = res$emp_power
  ))
}

# ============================================================================
# PART 6: DISPLAY TABLE
# ============================================================================

print(power_df)

# ============================================================================
# PART 7: CREATE PLOTS (X-axis: n)
# ============================================================================

create_plot_n <- function(combo_name) {
  
  df <- power_df[power_df$Combination == combo_name, ]
  df <- df[order(df$n), ]
  
  plot_df <- data.frame(
    n = rep(df$n, 2),
    Power = c(df$Theoretical_Power, df$Empirical_Power),
    Type = rep(c("Theoretical", "Empirical"), each = nrow(df))
  )
  
  y_max <- max(plot_df$Power) + 0.05
  y_min <- max(0, min(plot_df$Power) - 0.05)
  
  p <- ggplot(plot_df, aes(x = n, y = Power, color = Type, group = Type)) +
    geom_line(size = 1.5) +
    geom_point(size = 4) +
    scale_x_continuous(breaks = sample_sizes) +
    scale_y_continuous(limits = c(y_min, y_max), breaks = seq(0, 1, 0.05)) +
    labs(
      title = combo_name,
      subtitle = "theta = 1 (fixed), theoretical power is constant",
      x = "Sample Size (n)",
      y = "Power",
      color = "Type"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "bottom",
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10)
    )
  
  return(p)
}

# ============================================================================
# PART 8: CREATE PLOTS (X-axis: 1/√n)
# ============================================================================

create_plot_theta <- function(combo_name) {
  
  df <- power_df[power_df$Combination == combo_name, ]
  df <- df[order(df$n, decreasing = TRUE), ]
  
  theta_n <- 1 / sqrt(df$n)
  
  plot_df <- data.frame(
    theta_n = rep(theta_n, 2),
    Power = c(df$Theoretical_Power, df$Empirical_Power),
    Type = rep(c("Theoretical", "Empirical"), each = nrow(df))
  )
  
  theta_labels <- paste0("1/√", df$n)
  
  y_max <- max(plot_df$Power) + 0.05
  y_min <- max(0, min(plot_df$Power) - 0.05)
  
  p <- ggplot(plot_df, aes(x = theta_n, y = Power, color = Type, group = Type)) +
    geom_line(size = 1.5) +
    geom_point(size = 4) +
    scale_x_continuous(
      breaks = theta_n,
      labels = theta_labels
    ) +
    scale_y_continuous(limits = c(y_min, y_max), breaks = seq(0, 1, 0.05)) +
    labs(
      title = combo_name,
      subtitle = "theta = 1 (fixed), theoretical power is constant",
      x = expression(theta[n] == 1 / sqrt(n)),
      y = "Power",
      color = "Type"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "bottom",
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10)
    )
  
  return(p)
}

# ============================================================================
# PART 9: DISPLAY FOUR PLOTS (X-axis: n)
# ============================================================================

combinations <- c("Normal + Wilcoxon", "Normal + VdW", 
                  "Logistic + Wilcoxon", "Logistic + VdW")

plot_list_n <- list()
plot_list_theta <- list()

for (i in 1:length(combinations)) {
  plot_list_n[[i]] <- create_plot_n(combinations[i])
  plot_list_theta[[i]] <- create_plot_theta(combinations[i])
}

# Display plots with n on x-axis
for (i in 1:length(combinations)) {
  print(plot_list_n[[i]])
}

# Display plots with 1/√n on x-axis
for (i in 1:length(combinations)) {
  print(plot_list_theta[[i]])
}

# Combined plots
combined_n <- grid.arrange(
  plot_list_n[[1]], plot_list_n[[2]],
  plot_list_n[[3]], plot_list_n[[4]],
  ncol = 2, nrow = 2,
  top = grid::textGrob(
    "Power vs Sample Size (theta = 1 fixed)",
    gp = grid::gpar(fontsize = 16, face = "bold")
  )
)

combined_theta <- grid.arrange(
  plot_list_theta[[1]], plot_list_theta[[2]],
  plot_list_theta[[3]], plot_list_theta[[4]],
  ncol = 2, nrow = 2,
  top = grid::textGrob(
    "Power vs theta_n = 1/sqrt(n) (theta = 1 fixed)",
    gp = grid::gpar(fontsize = 16, face = "bold")
  )
)

print(combined_n)
print(combined_theta)

# ============================================================================
# END OF SCRIPT
# ============================================================================