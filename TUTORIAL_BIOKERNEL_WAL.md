<!--
# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

-->

# 🚀 BioKernel v3.1: The Workflow Abstraction Layer (WAL)

**Date:** February 14, 2026
**Architecture:** Platform-Agnostic Agentic OS

The **BioKernel** has evolved into a generalized **Workflow Abstraction Layer (WAL)**. 
It no longer depends on any single AI provider. Instead, it uses a plugin-based architecture to connect with **Gemini**, **OpenAI**, **Anthropic**, or **Local Models** (Llama 3 via Ollama) interchangeably.

## 🌟 Key Features

1.  **Provider Agnostic**: Switch between Gemini 2.0, GPT-4o, or Local LLMs just by changing a config file.
2.  **Standardized I/O**: All agents communicate using a unified `LLMRequest` / `LLMResponse` schema.
3.  **Resilient Factory**: The `LLMFactory` automatically loads the best available provider (e.g., falls back to Local if API keys are missing).
4.  **Extensible**: Add new models by implementing the `LLMProvider` interface.

## 🛠️ Configuration

Edit `platform/config.yaml` to choose your brain:

```yaml
# platform/config.yaml
provider:
  name: "gemini"  # Change to "local" or "openai"
```

## 🚀 How to Run

### 1. Set Credentials (Optional for Local)
```bash
export GOOGLE_API_KEY="your_key"
# OR
export OPENAI_API_KEY="your_key"
```

### 2. Start the Server
The BioKernel automatically detects the configured provider.
```bash
python3 platform/biokernel/server.py
```

*Output:*
- `🔌 [BioKernel] Loading Primary Provider: Gemini` (if key present)
- `🔌 [BioKernel] Loading Fallback Provider: Local` (if key missing)

### 3. Interact
The API remains the same, regardless of the underlying model.

```bash
curl -X POST "http://localhost:8000/v1/agent/run" \
     -H "Content-Type: application/json" \
     -d '{"query": "Analyze the role of TP53 in cancer."}'
```

## 📂 Architecture Map

- **`platform/interface/`**:
    - `llm_provider.py`: The Abstract Base Class (ABC) defining the contract.
- **`platform/schema/`**:
    - `io_types.py`: Standardized Pydantic models for Requests/Responses.
- **`platform/adapters/`**:
    - `factory.py`: The dynamic loader.
    - `gemini_adapter.py`: Google Gemini implementation.
    - `local_adapter.py`: Offline/Local model implementation.
- **`platform/biokernel/`**:
    - `server.py`: The OS Kernel using the WAL.

## 🔮 Adding a New Provider

1. Create `platform/adapters/my_adapter.py`.
2. Inherit from `LLMProvider`.
3. Implement `generate()` and `run_reasoning_loop()`.
4. Register it in `LLMFactory`.

<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->