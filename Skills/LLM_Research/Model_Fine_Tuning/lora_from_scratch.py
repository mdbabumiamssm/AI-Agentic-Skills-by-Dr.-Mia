import numpy as np

def lora_demo():
    """
    Demonstrates the mathematics of Low-Rank Adaptation (LoRA).
    
    Standard Fine-Tuning:
        W_new = W_old + delta_W
        (delta_W has same size as W_old, e.g., 1000x1000 = 1M params)
        
    LoRA:
        W_new = W_old + (B @ A)
        B is d x r
        A is r x k
        If r=4, params = 4000 + 4000 = 8k params (vs 1M!)
    """
    print("--- LoRA Mathematics Demo ---")
    
    # Dimensions
    d_in = 10
    d_out = 10
    rank = 2  # Low rank
    
    # 1. Pre-trained Weight Matrix (Frozen)
    W0 = np.random.randn(d_out, d_in)
    print(f"Pre-trained Weights (W0) shape: {W0.shape}")
    
    # 2. LoRA Matrices (Trainable)
    # A is usually initialized with Gaussian noise
    A = np.random.randn(rank, d_in)
    # B is initialized to ZEROS so training starts with no change
    B = np.zeros((d_out, rank))
    
    print(f"LoRA A shape: {A.shape}")
    print(f"LoRA B shape: {B.shape}")
    
    # 3. Input Vector
    x = np.random.randn(d_in, 1)
    
    # 4. Forward Pass
    # Standard path
    h_orig = W0 @ x
    
    # LoRA path: (B @ A) @ x
    # Note: We compute B@A implicitly or explicitly. 
    # Usually: B @ (A @ x) is faster.
    lora_delta = B @ (A @ x)
    
    h_final = h_orig + lora_delta
    
    print(f"\nOriginal Output (first 3): {h_orig.flatten()[:3]}")
    print(f"LoRA Delta (should be 0 initially): {lora_delta.flatten()[:3]}")
    print(f"Final Output: {h_final.flatten()[:3]}")
    
    # 5. Simulate Training (Update B slightly)
    print("\n--- Simulating 1 Step of Training ---")
    print("Updating Matrix B...")
    B += 0.1 # Arbitrary update
    
    lora_delta_new = B @ (A @ x)
    h_final_new = h_orig + lora_delta_new
    
    print(f"New LoRA Delta: {lora_delta_new.flatten()[:3]}")
    print(f"New Final Output: {h_final_new.flatten()[:3]}")
    print("The model output has changed without touching W0!")

if __name__ == "__main__":
    lora_demo()
