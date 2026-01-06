# LoRA (Low-Rank Adaptation)

**Category:** LLM Research / Fine-Tuning
**Difficulty:** Advanced

## Overview
Fine-tuning massive LLMs (70B+ parameters) is computationally expensive. **LoRA** freezes the pre-trained model weights and injects trainable rank decomposition matrices into each layer of the Transformer architecture.

## The Math
Instead of training the full weight update $\Delta W$, we approximate it with two smaller matrices $A$ and $B$:
$$ W_{new} = W_{0} + \Delta W = W_{0} + BA $$

Where:
-   $W_0 \in \mathbb{R}^{d \times k}$ (Frozen)
-   $B \in \mathbb{R}^{d \times r}$ (Trainable, initialized to 0)
-   $A \in \mathbb{R}^{r \times k}$ (Trainable, Gaussian init)
-   $r \ll \min(d, k)$ (Rank, e.g., 8 or 16)

## Benefits
1.  **Efficiency:** Reduces trainable parameters by up to 10,000x.
2.  **Storage:** Checkpoints are just a few MBs (only A and B).
3.  **Modularity:** You can swap "adapters" (e.g., a "Code Adapter" or a "Medical Adapter") at runtime without reloading the base model.

## Implementation (`lora_from_scratch.py`)
The script uses NumPy to demonstrate the exact matrix operations involved in the forward pass of a LoRA-adapted layer.

```