# Probability & Statistics for AI

**Category:** Mathematics
**Difficulty:** Beginner/Intermediate

## Overview
Machine Learning is essentially "Applied Statistics." Understanding distributions and hypothesis testing is crucial for:
1.  **Data Analysis:** Understanding the shape of your input features.
2.  **A/B Testing:** Validating if Model A is truly better than Model B.
3.  **Bayesian Methods:** Quantifying uncertainty in predictions.

## Key Concepts

### 1. Hypothesis Testing (T-Test)
Used to determine if there is a significant difference between the means of two groups.
-   **P-Value:** The probability of observing the data if the Null Hypothesis (no difference) were true.
-   **Threshold (Alpha):** Typically 0.05 (5%). If p < 0.05, we call the result "Statistically Significant."

### 2. Distributions
-   **Normal (Gaussian):** The "Bell Curve." Most natural phenomena follow this (CLT).
-   **Bernoulli/Binomial:** Coin flips (Success/Failure).
-   **Poisson:** Counts of rare events over time (e.g., server crashes, genetic mutations).

## Implementation (`stats_essentials.py`)
Demonstrates:
-   Generating synthetic data for a drug trial.
-   Running a T-Test using `scipy.stats`.
-   Simulating distributions using `numpy`.
