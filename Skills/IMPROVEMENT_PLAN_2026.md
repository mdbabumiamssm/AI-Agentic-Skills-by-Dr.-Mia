# Comprehensive Improvement Plan: Skills Repository (2026)

**Date:** January 5, 2026
**Status:** DRAFT
**Author:** AI Assistant

## 1. Executive Summary
The current `Skills` repository serves as a valuable **structural scaffold**, demonstrating the organization of AI agents for biomedical tasks. However, many implementations are currently **mocks**, "Lite" versions, or lack integration with production-grade libraries.

To align with the **2025/2026 State of the Art (SOTA)**, we must transition from "demonstration code" to "functional prototypes" that leverage:
1.  **Autonomous Agentic Workflows:** Moving beyond simple prompts to multi-step reasoning, planning, and tool use.
2.  **Real Scientific Libraries:** replacing mocks with `rdkit`, `biopython`, and `scikit-learn`.
3.  **Modern LLM Patterns:** Implementing Chain-of-Verification (CoV), MedPrompt, and Self-Refinement loops.

## 2. Architectural Upgrades (Computer Science & Agentic AI)

### A. Agent Framework Standardization
*   **Current:** Ad-hoc classes (`Agent`, `Supervisor`) with hardcoded "personas" and mock generation.
*   **Goal:** Create a robust `BaseAgent` class with:
    *   **LLM Adapter:** A unified interface for OpenAI, Anthropic, and local models (LLaMA/Mistral).
    *   **Tool Usage:** A standard decorator/protocol for defining tools agents can call.
    *   **Memory:** Vector-store integration (Milvus/Chroma) for long-term recall.

### B. Multi-Agent Orchestration
*   **Target:** `Agentic_AI/Multi_Agent_Systems/`
*   **Upgrade:** Implement a "Debate" and "Supervisor" pattern that uses **structured output** (JSON) to control flow, rather than simple string parsing.
*   **New Pattern:** "Language Agent Tree Search" (LATS) for complex problem solving.

## 3. Domain-Specific Enhancements

### A. Drug Discovery (`Drug_Discovery/`)
*   **ChemCrow Tools (`chem_tools.py`):**
    *   *Action:* Remove "Lite" mocks. Implement real `rdkit` functions for MolWt, LogP, TPSA, and QED.
    *   *Addition:* Add substructure matching (SMARTS) for toxicity alerts.
*   **AgentD:**
    *   *Action:* Implement a "Self-Driving Lab" interface mockup that simulates sending commands to a robotic platform (Opentrons).

### B. Clinical (`Clinical/`)
*   **Note Summarization:**
    *   *Upgrade:* Implement **MedPrompt** strategies (few-shot chain-of-thought with ensemble refinement).
    *   *Standard:* Structure outputs using FHIR-compliant JSON schemas.
*   **Trial Matching:**
    *   *Upgrade:* Add a "RAG" (Retrieval Augmented Generation) component to query a local database of `clinicaltrials.gov` JSONs.

## 4. Immediate Implementation Steps

1.  **Upgrade `debate_supervisor.py`:** Transform into a flexible Multi-Agent Orchestrator.
2.  **Upgrade `chem_tools.py`:** Inject real RDKit cheminformatics logic.
3.  **Create `medprompt_utils.py`:** Add a utility for advanced clinical prompting.
