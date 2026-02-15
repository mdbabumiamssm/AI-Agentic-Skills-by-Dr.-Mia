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

# Distributed Systems for AI

**Category:** Computer Science / Systems
**Difficulty:** Advanced

## Overview
Distributed systems are the backbone of modern AI. Training Large Language Models (LLMs) requires thousands of GPUs working in parallel. This module demonstrates the **Parameter Server** and **Actor** patterns, which are central to frameworks like **Ray** and **PyTorch Distributed**.

## Key Concepts

### 1. The CAP Theorem
In a distributed data store, you can only provide two of the three guarantees:
-   **Consistency:** Every read receives the most recent write or an error.
-   **Availability:** Every request receives a (non-error) response, without the guarantee that it contains the most recent write.
-   **Partition Tolerance:** The system continues to operate despite an arbitrary number of messages being dropped (or delayed) by the network between nodes.

### 2. Parallelism Strategies
*   **Data Parallelism:** The model is replicated on every device. Each device processes a different chunk of data (batch). Gradients are aggregated (averaged) to update the model.
*   **Model Parallelism:** The model is too large to fit on one device. Different layers (or parts of layers) effectively live on different devices.

### 3. The Actor Model (Ray)
Instead of just sending passive data, we send active "Actors" (stateful worker processes) to different nodes. 
-   **Stateless Tasks:** Functions that take input and return output (e.g., data preprocessing).
-   **Stateful Actors:** Classes that maintain state (e.g., a Neural Network holding weights).

## Implementation (`ray_mock.py`)
This script simulates a `Parameter Server` architecture using a mock version of the Ray API.
1.  **Workers:** Calculate gradients on their local data.
2.  **Server:** Aggregates gradients and updates global weights.
3.  **Futures:** The `.remote()` call returns a future object, allowing asynchronous execution (simulated).


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->