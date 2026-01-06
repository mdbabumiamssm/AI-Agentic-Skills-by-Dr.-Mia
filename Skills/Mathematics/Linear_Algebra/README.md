# Linear Algebra for AI

**Category:** Mathematics
**Difficulty:** Beginner/Intermediate

## Overview
Linear Algebra is the language of machine learning. Data is represented as vectors and matrices, and operations on them are fundamental to neural networks.

## Key Concepts

### 1. Vectors and Matrices
- **Scalar:** A single number (0-tensor).
- **Vector:** An array of numbers (1-tensor). Represents a point in space or a direction.
- **Matrix:** A 2D array of numbers (2-tensor). Represents a linear transformation.

### 2. Matrix Multiplication (Dot Product)
Unlike element-wise multiplication, the dot product combines rows of the first matrix with columns of the second.
If $A$ is $(m \times n)$ and $B$ is $(n \times p)$, then $C = AB$ is $(m \times p)$.

### 3. Eigenvalues and Eigenvectors
For a square matrix $A$, an eigenvector $v$ is a vector that does not change direction when $A$ is applied to it. It only scales by a factor $\lambda$ (the eigenvalue).
$$ Av = \lambda v $$

## Usage
This module requires `numpy`.

```bash
pip install numpy
python matrix_ops.py
```

```