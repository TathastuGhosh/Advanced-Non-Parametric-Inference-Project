<div align="center">

# Application of Le Cam's Lemma in Two-Sample Location Problem and Its Finite-Sample Performance

### Master's Dissertation

**Tathastu Ghosh**

*M.Sc. Statistics*  
**Department of Statistics**  
**University of Calcutta**

<br>

[![R](https://img.shields.io/badge/R-4.0+-276DC3?style=for-the-badge&logo=r)](https://www.r-project.org/)
![Master's Dissertation](https://img.shields.io/badge/Master's-Dissertation-8A2BE2?style=for-the-badge)
![Monte Carlo Simulation](https://img.shields.io/badge/Monte%20Carlo-Simulation-009688?style=for-the-badge)
![Nonparametric Statistics](https://img.shields.io/badge/Nonparametric-Statistics-E91E63?style=for-the-badge)
![University of Calcutta](https://img.shields.io/badge/University%20of-Calcutta-7B1FA2?style=for-the-badge)

</div>

---

## 📖 Overview

This repository contains the complete R code, simulation scripts, and visualization tools for the Master's dissertation titled **"Application of Le Cam's Lemma in Two-Sample Location Problem and Its Finite-Sample Performance"** . The project investigates the asymptotic behavior of linear rank test statistics under contiguous local alternatives through extensive Monte Carlo simulations.

The study validates Le Cam's lemma by examining:
- **Distributional Validation**: Assessing how well the asymptotic normal distribution approximates finite-sample distributions under local alternatives
- **Power Analysis (First Perspective)**: Evaluating the accuracy of asymptotic power approximations by varying the effect size \(\theta\) for fixed sample sizes
- **Power Analysis (Second Perspective)**: Examining the convergence of empirical power to the asymptotic power by fixing the effect size θ = 1 and varying the sample size n

### Key Features

- ✅ Comprehensive simulation framework for two-sample rank tests
- ✅ Implementation of Wilcoxon and Van der Waerden scores
- ✅ Support for Normal and Logistic distributions
- ✅ Kolmogorov-Smirnov distance computation for distributional validation
- ✅ Theoretical and empirical power curve generation (two perspectives)
- ✅ Publication-ready plots and tables
- ✅ Fully reproducible results with fixed random seeds

---

## 📊 Simulation Design

### Parameters

| Parameter | Levels | Description |
|:---|:---|:---|
| **Sample Size (n)** | 30, 50, 70, 100, 200 | Total sample size (balanced design, n₁ = n₂ = n/2) |
| **Distributions** | Normal, Logistic | Underlying population distributions |
| **Score Functions** | Wilcoxon, Van der Waerden | Score functions for rank tests |
| **Local Shift** | θ = 1 | Fixed local alternative parameter |
| **Replications** | 10,000 (Distributional) / 5,000 (Power) | Monte Carlo replications per scenario |

### Data Generation

Under the local alternative θₙ = θ/√n

- **Normal Distribution**:
  - Group 1: X ∼ N(0, 1)
  - Group 2: Y ∼ N(θₙ, 1)

- **Logistic Distribution**:
  - Group 1: X ∼ Logistic(0, 1)
  - Group 2: Y ∼ Logistic(θₙ, 1)

### Key Theoretical Parameters

| Distribution | Score | g₀ | σ²(ψ) |
|:---|:---|:---|:---|
| Normal | Wilcoxon | 0.4886 | 1 |
| Normal | **Van der Waerden** | **1.0000** | **1** | 
| Logistic | **Wilcoxon** | **0.2887** | **1/3** | 
| Logistic | Van der Waerden | 0.5642 | 1/3 | 

---

## 📈 Results Summary

### Distributional Validation (KS Distances)

The Kolmogorov-Smirnov distances confirm the theoretical optimality of score functions:

| Distribution | Score | n=30 | n=50 | n=70 | n=100 | n=200 |
|:---|:---|:---|:---|:---|:---|:---|
| Normal | Wilcoxon | 0.0475 | 0.0358 | 0.0282 | 0.0215 | 0.0118 |
| Normal | **Van der Waerden** | **0.0324** | **0.0241** | **0.0186** | **0.0132** | **0.0075** |
| Logistic | **Wilcoxon** | **0.1183** | **0.0952** | **0.0825** | **0.0684** | **0.0452** |
| Logistic | Van der Waerden | 0.1856 | 0.1687 | 0.1554 | 0.1432 | 0.1258 |

### Power Analysis - First Perspective (Fixed n, Varying θ)

The hierarchy of power is consistently maintained across all sample sizes:


Normal + VdW > Logistic + Wilcoxon > Normal + Wilcoxon > Logistic + VdW


### Power Analysis - Second Perspective (Fixed θ = 1, Varying n)

The convergence of empirical power to the asymptotic power reveals:

| Combination | Convergence Rate |
|:------------|:-----------------|
| **Normal + Van der Waerden** | **Fastest** |
| Logistic + Wilcoxon | Fast |
| Normal + Wilcoxon | Slow |
| **Logistic + Van der Waerden** | **Slowest** |

---

## Key Findings

1. **Optimality Confirmed**: The Van der Waerden score is optimal for Normal data; the Wilcoxon score is optimal for Logistic data.

2. **Misspecification Penalty**: Using a sub-optimal score function results in substantial power loss (40–50% for large effect sizes).

3. **Convergence**: Empirical distributions approach asymptotic normal approximations as sample size increases.

4. **Heavier Tails**: Logistic distribution exhibits slower convergence and lower power ceilings compared to Normal distribution.

5. **Rate of Convergence**: The rate at which empirical power approaches the asymptotic limit depends critically on the alignment between the score function and the underlying distribution. Optimal combinations converge rapidly, while sub-optimal combinations exhibit persistent gaps even at larger sample sizes.

---

## 🛠️ Required R Packages

```r
install.packages(c(
    "ggplot2",      # Data visualization
    "dplyr",        # Data manipulation
    "tidyr",        # Data reshaping
    "gridExtra",    # Arranging multiple plots
    "knitr"         # Table formatting
))