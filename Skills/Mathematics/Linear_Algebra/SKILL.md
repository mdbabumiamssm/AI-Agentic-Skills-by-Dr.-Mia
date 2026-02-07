---
name: 'tensor-operations'
description: 'Tensor Operations'
---


# Tensor Operations

Fundamental linear algebra operations for understanding Transformers and attention mechanisms.

## When to Use This Skill

*   When you need to manually compute Attention mechanisms.
*   For educational purposes to understand "Self-Attention".
*   To creating custom masks for sequence modeling.

## Core Capabilities

1.  **Scaled Dot Product Attention**: `softmax(QK^T / sqrt(d_k))`.
2.  **Causal Masking**: Create lower-triangular masks for GPT-style generation.

## Workflow

1.  **Input**: Query, Key, Value matrices.
2.  **Execute**: Run the script.
3.  **Output**: Attention scores and weighted values.

## Example Usage

**User**: "Compute attention for these matrices."

**Agent Action**:
```bash
python3 Skills/Mathematics/Linear_Algebra/tensor_operations.py
```
