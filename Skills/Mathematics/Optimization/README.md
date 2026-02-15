<!--
# COPYRIGHT NOTICE
# This file is part of the "Universal AI Agentic Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

-->

# Optimization: Gradient Descent

**Category:** Mathematics / Machine Learning
**Difficulty:** Beginner

## Overview
Gradient Descent is the primary optimization algorithm used to train machine learning models and neural networks. It iteratively minimizes a function by moving in the direction of the steepest descent (the negative of the gradient).

## The Math
For a function $J(\theta)$ (the Cost Function), the update rule for parameter $\theta$ is:

$$ \theta_{new} = \theta_{old} - \alpha \cdot \nabla J(\theta) $$

Where:
*   $\alpha$ (alpha) is the **Learning Rate** (step size).
*   $\nabla J(\theta)$ is the **Gradient** (derivative) of the cost function at $\theta$.

## Code Example
The script `gradient_descent.py` minimizes the quadratic function:
$$ f(x) = x^2 + 4x + 4 $$

We know from calculus that the minimum is at $x = -2$. The script demonstrates how the algorithm starts at a random point (e.g., $x=5$) and "walks" down the curve to find -2.

## Key Concepts
1.  **Learning Rate:** If too high, it overshoots. If too low, it's too slow.
2.  **Epochs:** The number of iteration steps.
3.  **Convexity:** This simple example uses a convex function (bowl shape), guaranteeing a global minimum. Neural networks involve non-convex functions.



<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->