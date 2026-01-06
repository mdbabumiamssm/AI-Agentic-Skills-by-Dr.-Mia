# Memory Systems for Agents

**Category:** Agentic AI / Architecture
**Difficulty:** Intermediate

## Overview
A stateless agent typically "forgets" everything after the window closes. A robust Memory System allows agents to maintain continuity, personalize interactions, and learn over time.

## Architecture

### 1. Short-Term Memory (Working Memory)
-   **Content:** The current conversation context window.
-   **Mechanism:** Direct prompt injection.
-   **Limit:** Constrained by the model's context length (e.g., 128k tokens).

### 2. Episodic Memory (Experience)
-   **Content:** Past interactions, events, and user history. "I remember we discussed Python last Tuesday."
-   **Mechanism:** Time-series logs, often summarized periodically.
-   **Storage:** JSON, SQL, or Vector DB (for similarity search).

### 3. Semantic Memory (Knowledge)
-   **Content:** General facts about the world or the user. "Alice is a Python developer."
-   **Mechanism:** Embeddings stored in a Vector Database.
-   **Process:** RAG (Retrieval Augmented Generation).

## Implementation (`memory_architecture.py`)
This script demonstrates a dual-memory system:
-   **Episodic:** Logs every turn with a timestamp.
-   **Semantic:** "Learns" facts (heuristically) and stores them in a mock Vector Store for recall.
