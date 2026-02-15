# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import numpy as np

def basic_operations():
    print("--- Basic Operations ---")
    A = np.array([[1, 2], [3, 4]])
    B = np.array([[5, 6], [7, 8]])
    
    print(f"Matrix A:\n{A}")
    print(f"Matrix B:\n{B}")
    
    # Element-wise addition
    print(f"A + B:\n{A + B}")
    
    # Element-wise multiplication
    print(f"A * B (Hadamard Product):\n{A * B}")
    
    # Matrix Multiplication (Dot Product)
    print(f"A @ B (Dot Product):\n{A @ B}")

def transformations():
    print("\n--- Linear Transformations ---")
    v = np.array([1, 0]) # Unit vector on X axis
    print(f"Vector v: {v}")
    
    # Rotation Matrix (90 degrees counter-clockwise)
    theta = np.radians(90)
    R = np.array([
        [np.cos(theta), -np.sin(theta)],
        [np.sin(theta), np.cos(theta)]
    ])
    
    v_rotated = R @ v
    print(f"Rotation Matrix (90 deg):\n{R}")
    print(f"Rotated v: {v_rotated}")

def eigen_decomposition():
    print("\n--- Eigenvalues and Eigenvectors ---")
    # A symmetric matrix ensures real eigenvalues
    M = np.array([[4, 1], [1, 2]]) 
    print(f"Matrix M:\n{M}")
    
    vals, vecs = np.linalg.eig(M)
    
    for i in range(len(vals)):
        print(f"Eigenvalue {i}: {vals[i]:.2f}")
        print(f"Eigenvector {i}: {vecs[:, i]}")
        
    print("Recall: M @ v = lambda * v")
    v0 = vecs[:, 0]
    lambda0 = vals[0]
    print(f"M @ v0: {M @ v0}")
    print(f"lambda0 * v0: {lambda0 * v0}")

if __name__ == "__main__":
    basic_operations()
    transformations()
    eigen_decomposition()

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
