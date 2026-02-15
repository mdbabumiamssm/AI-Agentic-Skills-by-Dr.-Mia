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

# Universal AI Agentic Skills (2026)

![Status](https://img.shields.io/badge/Status-Active-green)
![Architecture](https://img.shields.io/badge/Architecture-Universal%20Agentic%20OS-blueviolet)
![Domain](https://img.shields.io/badge/Domain-Universal%20%7C%20Finance%20%7C%20Coding%20%7C%20Science-orange)
![Tech](https://img.shields.io/badge/Tech-Gemini%202.0%20%7C%20OpenAI%20o3%20%7C%20Claude%203.7%20%7C%20MCP-blue)

> **⚠️ IMPORTANT DISCLAIMER & COPYRIGHT NOTICE**
> 
> This repository, its architecture, agent designs, and specific implementations are the intellectual property of **MD BABU MIA, PhD**.
> 
> While open-source components are licensed under MIT, the unique curation, "Biomedical OS" architecture, and agentic workflows are proprietary to the author. 
> 
> **If you fork, clone, or copy this repository for public use, you MUST:**
> 1.  Retain this copyright notice.
> 2.  Explicitly credit **MD BABU MIA, PhD** as the original author.
> 3.  Link back to the original repository. 
> 
> *Plagiarism or uncredited redistribution is strictly prohibited.*

---

## 🚀 The Universal Operating System for AI Agents

**One Platform. Infinite Agents. Zero Limits.**

This is not just a collection of scripts; it is a **Universal Agentic Operating System (OS)** designed to empower professionals across every high-demand industry. By abstracting the "brain" (LLM) from the "body" (Tools), we allow you to deploy state-of-the-art autonomous agents for **Coding, Finance, Research, and Science** instantly.

Whether you need to refactor a 100k-line codebase, predict market volatility, simulate quantum docking, or design a clinical trial, this platform has a specialized agent ready for you.

---

## 🎯 Who is this for? (2025-2026 Industry Trends)

We have organized 600+ skills into specialized domains targeting the most critical needs of the modern workforce.

### 💻 For Software Engineers & Architects
*   **Automated Refactoring:** Deploy `Codebase Investigator` to map legacy repositories, identify debt, and generate modernization plans automatically.
*   **Spec-Driven Development:** Agents that convert high-level PRDs (Product Requirement Docs) into executable boilerplate code.
*   **DevOps & CI/CD:** Autonomous agents that monitor pipelines, auto-fix linting errors, and generate infrastructure-as-code (Terraform/Pulumi).
*   **Key Skills:** `Software_Engineering`, `Computer_Science`, `Workflow_Management`.

### 🔬 For Scientists & Physicists
*   **Quantum Biotech:** Bridge the gap between physics and biology with agents for molecular dynamics simulations and quantum docking.
*   **Lab Automation:** "Self-Driving Lab" agents that design experiments, control robotic liquid handlers (Opentrons), and analyze results loops.
*   **Mathematical Proofs:** Specialized agents for symbolic math, differential equations, and complex optimization problems.
*   **Key Skills:** `Quantum_Biotech`, `Mathematics`, `Structural_Biology`, `Lab_Automation`.

### 🤖 For AI Researchers
*   **Swarm Architecture:** Ready-to-use templates for "Plan-and-Solve," "ReAct," and "Map-Reduce" agent topologies.
*   **Model Evaluation:** Automated benchmarking pipelines to test LLMs against domain-specific datasets (Medical, Coding, Math).
*   **Prompt Optimization:** The `Meta-Prompter` module that uses AI to rewrite and optimize its own instructions for peak performance.
*   **Key Skills:** `Agentic_AI`, `LLM_Research`, `Foundation_Models`.

### 🧬 For Biomedical Experts (The Core)
*   **Genomics:** End-to-end pipelines for Single-Cell RNA-seq, Spatial Transcriptomics, and Variant Calling.
*   **Clinical:** Decision support systems, automated prior authorization, and "ChatEHR" for patient data interaction.
*   **Pharma:** Generative antibody design (MAGE), small molecule evolution, and regulatory submission drafting.
*   **Key Skills:** `Genomics`, `Clinical`, `Drug_Discovery`.

---

## 🌟 Architecture: The "Bio-OS" Kernel

At the heart of this system is the **BioKernel**, a high-performance runtime environment.

*   **Workflow Abstraction Layer (WAL):** A plugin-based system that allows you to swap the "brain" of the agent. 
    *   *Need Speed?* Use **Gemini 2.0 Flash**.
    *   *Need Reasoning?* Use **OpenAI o3** or **Claude 3.7**.
    *   *Need Privacy?* Use **Local Llama 3**.
*   **Universal Skill Definition (USDL):** All skills are defined in a standardized JSON/YAML schema, making them portable across any LLM framework (LangChain, Semantic Kernel, AutoGen).

### 💰 For Finance & Quant Professionals
*   **Market Swarms:** Multi-agent systems that aggregate news, sentiment, and technical indicators to signal market shifts in real-time.
*   **Automated Compliance:** Agents that audit transactions against regulatory frameworks (SOX, GDPR) and generate risk reports.
*   **Portfolio Management:** `Finance_Agent_Team` capable of rebalancing portfolios and simulating stress tests based on macro-economic scenarios.
*   **Key Skills:** `Data_Science`, `Mathematics`, `External_Collections/ai_finance_agent`

---

## 📂 Global Capability Map

```text
Skills/
├── Agentic_AI/           # 🧠 The Brains: Swarms, Planners, Memory Systems
├── Computer_Science/     # 💻 The Engineering: Algos, Distributed Systems
├── Data_Science/         # 📊 The Analytics: Visualization, ETL Pipelines
├── Finance_Agents/       # 💰 The Markets: Trading, Risk, Compliance (External)
├── Mathematics/          # 🧮 The Logic: Optimization, Linear Algebra
├── Quantum_Biotech/      # ⚛️ The Physics: Simulation, Docking
├── Software_Engineering/ # 🛠️ The Tools: Refactoring, Testing, DevOps
├── Clinical/             # 🏥 The Hospital: EHR, Diagnosis, Trials
├── Genomics/             # 🧬 The Lab: Sequencing, CRISPR, Single-Cell
└── Drug_Discovery/       # 💊 The Cure: Chemistry, Antibodies, Pharma
```

## 🛠️ Quick Start

### 1. Select Your Brain
Edit `platform/config.yaml` to choose your intelligence provider.

```yaml
provider:
  name: "gemini" # Options: gemini, openai, anthropic, local
```

### 2. Boot the OS
Initialize the kernel to load the 600+ skills into active memory.

```bash
python3 platform/biokernel/server.py
```

### 3. Command Your Army
Send natural language commands to the Universal Agent.

**"Analyze this repository and refactor the auth module:"**
```bash
curl -X POST "http://localhost:8000/v1/agent/run" \
     -d '{"query": "Analyze ./src/auth and propose a refactor pattern."}'
```

**"Find an arbitrage opportunity in these datasets:"**
```bash
curl -X POST "http://localhost:8000/v1/agent/run" \
     -d '{"query": "Compare dataset A and B for pricing anomalies."}'
```

**"Design a guide RNA for this gene target:"**
```bash
curl -X POST "http://localhost:8000/v1/agent/run" \
     -d '{"query": "Design CRISPR gRNA for TP53 exon 4."}'
```

---

## 👤 Author & Maintainer

**MD BABU MIA, PhD**  
*Assistant Professor of Hematology & Medical Oncology, Machine Learning -AI*  
Mount Sinai Tisch Cancer Institute  
Icahn School of Medicine at Mount Sinai  
New York, NY 10029  
**Email:** md.babu.mia@mssm.edu  

---

## 📄 License

**Copyright (c) 2026 MD BABU MIA, PhD.**  
All rights reserved.

This project is licensed under the MIT License for open-source components, but the unique architectural design and agentic workflows are the intellectual property of the author. **Attribution is mandatory.**

<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->
