# =============================================================================
# APPLICATION OF LE CAM'S LEMMA: FINITE-SAMPLE PERFORMANCE
# =============================================================================
# This script simulates linear rank statistics under local alternatives
# and compares empirical distributions with theoretical normal approximations.
# =============================================================================

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

# Set seed for reproducibility
set.seed(123)

# =============================================================================
# 1. HELPER FUNCTIONS
# =============================================================================

# Compute Wilcoxon scores
wilcoxon_scores <- function(n) {
  return(1:n)
}

# Compute Van der Waerden scores
vdw_scores <- function(n) {
  return(qnorm((1:n) / (n + 1)))
}

# Compute linear rank statistic T_n
compute_Tn <- function(x, y, score_func) {
  n1 <- length(x)
  n2 <- length(y)
  n <- n1 + n2
  
  pooled <- c(x, y)
  ranks <- rank(pooled)
  c_i <- c(rep(0, n1), rep(1, n2))
  c_bar <- mean(c_i)
  a_ranks <- score_func(n)[ranks]
  T_n <- sum((c_i - c_bar) * a_ranks)
  
  return(T_n)
}

# Generate data under alternative (theta_n = 1/sqrt(n))
generate_alternative_data <- function(n1, n2, dist_type) {
  n <- n1 + n2
  theta_n <- 1 / sqrt(n)
  
  if (dist_type == "normal") {
    x <- rnorm(n1, mean = 0, sd = 1)
    y <- rnorm(n2, mean = theta_n, sd = 1)
  } else if (dist_type == "logistic") {
    x <- rlogis(n1, location = 0, scale = 1)
    y <- rlogis(n2, location = theta_n, scale = 1)
  }
  
  return(list(x = x, y = y))
}

# Compute theoretical parameters under alternative
theoretical_params <- function(n1, n2, dist_type, score_type) {
  n <- n1 + n2
  
  if (dist_type == "normal") {
    sigma2_psi <- 1
    if (score_type == "wilcoxon") {
      g0 <- 1 / (2 * sqrt(pi))
    } else if (score_type == "vdw") {
      g0 <- 1
    }
  } else if (dist_type == "logistic") {
    sigma2_psi <- 1/3
    if (score_type == "wilcoxon") {
      g0 <- 1/6
    } else if (score_type == "vdw") {
      g0 <- 1/sqrt(pi)
    }
  }
  
  mean_T <- (n1 * n2 / n) * g0 / sqrt(n)
  var_T <- (n1 * n2 / n) * sigma2_psi
  sd_T <- sqrt(var_T)
  
  return(list(mean = mean_T, sd = sd_T))
}

# =============================================================================
# 2. MAIN SIMULATION FUNCTION
# =============================================================================

run_simulation <- function(n, n_sim = 2000) {
  n1 <- n2 <- n/2
  
  cat("\nRunning simulation for n =", n, "\n")
  
  # Results storage
  all_results <- data.frame()
  
  # Score function definitions
  score_functions <- list(
    "Wilcoxon" = wilcoxon_scores,
    "VdW" = vdw_scores
  )
  
  score_types <- list(
    "Wilcoxon" = "wilcoxon",
    "VdW" = "vdw"
  )
  
  # Loop over score types
  for (score_name in names(score_functions)) {
    score_func <- score_functions[[score_name]]
    score_type <- score_types[[score_name]]
    
    # Loop over data distributions
    for (dist_type in c("Normal", "Logistic")) {
      dist_type_lower <- tolower(dist_type)
      
      # Storage for T_n values
      Tn_values <- numeric(n_sim)
      
      # Generate data and compute T_n
      for (i in 1:n_sim) {
        data <- generate_alternative_data(n1, n2, dist_type_lower)
        Tn_values[i] <- compute_Tn(data$x, data$y, score_func)
      }
      
      # Compute theoretical parameters
      theo_params <- theoretical_params(n1, n2, dist_type_lower, score_type)
      
      # Store results
      temp_results <- data.frame(
        Tn = Tn_values,
        n = rep(as.numeric(n), n_sim),
        score = rep(as.character(score_name), n_sim),
        distribution = rep(as.character(dist_type), n_sim),
        theo_mean = rep(as.numeric(theo_params$mean), n_sim),
        theo_sd = rep(as.numeric(theo_params$sd), n_sim)
      )
      all_results <- rbind(all_results, temp_results)
    }
  }
  
  # =========================================================================
  # 3. CREATE DENSITY PLOTS FOR THIS SAMPLE SIZE
  # =========================================================================
  
  plot_list <- list()
  
  for (score_name in c("Wilcoxon", "VdW")) {
    
    # Subset data for this score
    score_data <- subset(all_results, score == score_name)
    
    # Get theoretical parameters
    theo_mean <- unique(score_data$theo_mean)[1]
    theo_sd <- unique(score_data$theo_sd)[1]
    
    # Subset data by distribution
    normal_data <- subset(score_data, distribution == "Normal")
    logistic_data <- subset(score_data, distribution == "Logistic")
    
    # Create the plot
    p <- ggplot(score_data, aes(x = Tn)) +
      # Theoretical normal density
      stat_function(
        fun = function(x) dnorm(x, mean = theo_mean, sd = theo_sd),
        aes(color = "Theoretical"),
        linewidth = 1.2,
        linetype = "solid"
      ) +
      # Empirical density for Normal data
      geom_density(
        data = normal_data,
        aes(color = "Normal Data", fill = "Normal Data"),
        alpha = 0.3,
        linewidth = 0.8
      ) +
      # Empirical density for Logistic data
      geom_density(
        data = logistic_data,
        aes(color = "Logistic Data", fill = "Logistic Data"),
        alpha = 0.3,
        linewidth = 0.8
      ) +
      # Labels and theme
      labs(
        title = paste(score_name, "Scores"),
        x = expression(T[n]),
        y = "Density",
        color = "",
        fill = ""
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        legend.position = "bottom",
        legend.box = "horizontal",
        panel.grid.minor = element_blank()
      ) +
      scale_color_manual(
        values = c("Theoretical" = "black", 
                   "Normal Data" = "#1f78b4", 
                   "Logistic Data" = "#e31a1c")
      ) +
      scale_fill_manual(
        values = c("Normal Data" = "#1f78b4", 
                   "Logistic Data" = "#e31a1c")
      )
    
    plot_list[[score_name]] <- p
  }
  
  # Combine plots for this sample size
  combined_plot <- grid.arrange(
    plot_list[["Wilcoxon"]], 
    plot_list[["VdW"]],
    ncol = 2,
    top = grid::textGrob(
      paste("n =", n),
      gp = grid::gpar(fontsize = 14, fontface = "bold")
    )
  )
  
  return(list(
    combined_plot = combined_plot,
    results = all_results
  ))
}

# =============================================================================
# 4. RUN SIMULATIONS FOR ALL SAMPLE SIZES
# =============================================================================

sample_sizes <- c(30, 50, 70, 100, 200)

# Store all results
all_sim_results <- list()

for (n in sample_sizes) {
  # Adjust simulation iterations based on sample size
  n_sim <- ifelse(n <= 50, 3000, ifelse(n <= 100, 2000, 1000))
  
  result <- run_simulation(n, n_sim = n_sim)
  
  all_sim_results[[as.character(n)]] <- result
  
  # Display the plot
  print(result$combined_plot)
}
