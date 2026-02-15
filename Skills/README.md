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

# Skills Repository (2026 Edition)

> **The "Microservices Moment" for Biomedical AI Agents.**

![Status](https://img.shields.io/badge/Status-Active-green)
![Agents](https://img.shields.io/badge/Agents-Orchestrated-blue)
![Domain](https://img.shields.io/badge/Domain-Biotech%20%7C%20Clinical%20%7C%20Genomics-purple)
![Tech](https://img.shields.io/badge/Tech-MCP%20%7C%20DeepSeek%20%7C%20Gemini-orange)

## 🚀 Overview

This repository is a comprehensive library of **skills, agents, and mathematical foundations** for modern (2026) Artificial Intelligence. Unlike standard chatbot repos, this project focuses on **Agentic Workflows**—where autonomous systems plan, execute, use tools, and correct themselves to solve complex scientific problems.

We have aligned this codebase with the **State of the Art (SOTA) for 2026**, integrating Agentic patterns, Model Context Protocol (MCP), and rigorous scientific simulations.

## 🌟 Key Capabilities (New for 2026)

### 🧬 Genomics & Single Cell (New!)
*   **Universal Annotator:** `Genomics/Single_Cell/Cell_Type_Annotation/RNA/universal_annotator.py` wraps Marker-based, Deep Learning (CellTypist), and LLM-based annotation strategies.
*   **Pathway Scoring:** `Genomics/Single_Cell/Pathway_Analysis/sc_pathway_scorer.py` implements AUCell-like scoring for functional enrichment.
*   **Cell-Cell Comms:** `Genomics/Single_Cell/Cell_Cell_Communication/interaction_inference.py` infers Ligand-Receptor networks.
*   **Database:** A curated [Tool Database](Genomics/Single_Cell/Tool_Database.md) of 2026 single-cell tools (MultiKano, scPS, etc.).

### 🧠 Agentic AI (The Brain)
*   **Self-Correction:** `Agentic_AI/Agent_Architectures/Self_Correction/self_correction_agent.py` implements a Reflexion pattern for iterative improvement.
*   **Orchestrated Swarms:** `Agentic_AI/Multi_Agent_Systems/orchestrator.py` implements a Supervisor pattern that delegates tasks to specialized sub-agents.
*   **Plan-and-Solve:** `Agentic_AI/Agent_Architectures/Plan_and_Solve/` breaks down complex user queries into Directed Acyclic Graphs (DAGs).
*   **Async Runtime:** `Computer_Science/Distributed_Systems/agent_concurrency.py` provides a Ray-like async runtime for parallel agents.

### 🔌 Model Context Protocol (MCP) & Platform
*   **BioMCP Server:** `MCP_Servers/BioMCP/bio_mcp_server.py` implements a compliant MCP server exposing bio-tools.
*   **BioKernel Engine:** `../platform/biokernel/workflow_engine.py` orchestrates enterprise-grade agent workflows (Mining -> Design -> Safety) via FastAPI.
*   **Runtime Adapter:** `../platform/adapters/runtime_adapter.py` provides a unified execution layer, supporting intelligent mock simulations and real API calls.

### 🏥 Clinical & Operations (New!)
*   **Trial Matching:** `Clinical/Trial_Matching/trial_matching_agent.py` matches patient profiles to clinical trials using intelligent criteria mapping (LLM-driven).
*   **Prior Auth Appeals:** `Clinical/Prior_Authorization/appeals_agent.py` uses self-correction to iteratively refine arguments for overturning insurance denials.
*   **EHR/FHIR Integration:** `Clinical/EHR_FHIR_Integration/fhir_client.py` provides tools to search and retrieve patient data from FHIR R4 servers.
*   **Clinical NLP:** `Clinical/Clinical_NLP/entity_extractor.py` extracts medical entities (Diseases, Meds) from unstructured text.
*   **Opentrons Agent:** `Lab_Automation/Opentrons_Agent/opentrons_generator.py` generates liquid handling protocols from high-level intent.

### 💊 Drug Discovery & Genomics (Updated)
*   **Literature Mining:** `Research_Tools/Literature_Mining/mining_agent.py` uses the Runtime Adapter to extract novel targets from (simulated) texts.
*   **Molecule Evolution:** `Drug_Discovery/Molecule_Design/evolution_agent.py` designs de novo drugs with 'medicinal chemist' feedback loop.
*   **Safety Officer:** `Clinical/Safety/safety_agent.py` audits outputs for compliance and toxicity using semantic analysis.
*   **Variant Interpretation:** `Genomics/Variant_Interpretation/acmg_classifier.py` classifies genetic variants and generates AI-powered clinical reports.
*   **ChemCrow Tools:** `Drug_Discovery/ChemCrow_Tools/chem_tools.py` enables agents to calculate molecular properties (LogP, TPSA) and screen for toxicity.
*   **CRISPR Design:** `Genomics/CRISPR_Design_Agent/crispr_designer.py` automates gRNA selection and efficiency scoring for gene editing.
*   **Protein Structure:** `Drug_Discovery/Protein_Structure/esmfold_client.py` mocks ESMFold/AF3 inference for 3D structure prediction.

### 🧪 Clinical Simulators & Research
*   **Adaptive Clinical Trials:** `Clinical/Clinical_Trials/Adaptive_Trial_Design_Agent/adaptive_trial_sim.py` runs Bayesian MAMS simulations.
*   **MedPrompt:** `LLM_Research/Prompt_Engineering/medprompt.py` implements Microsoft's SOTA clinical reasoning strategy.
*   **Self-Driving Labs:** `Mathematics/Probability_Statistics/bayesian_optimization.py` enables autonomous experiment selection using Gaussian Processes.

### 📐 Math & CS (The Foundation)
*   **Tensor Operations:** `Mathematics/Linear_Algebra/tensor_operations.py` breaks down the math behind Attention mechanisms.
*   **Graph RAG:** `Computer_Science/Graph_Algorithms/knowledge_graph.py` provides traversal for Drug-Target-Disease interactions.

### 💻 Software Engineering (New!)
*   **React & Next.js Best Practices:** `Software_Engineering/Web_Development/` contains standardized rules (`SKILL.md`) for building modern, performant web UIs.
*   **Data Science Standards:** `Software_Engineering/Data_Science/Python_Pandas_Best_Practices/` provides guidelines for vectorized, memory-efficient data manipulation.
*   **Core Python:** `Software_Engineering/Core_Python_Best_Practices/` enforces modern typing and idiomatic Python 3.10+ patterns.

## 🤝 Dual Health Stacks (New)

### OpenAI Health Stack
*   **Care Copilot:** `Consumer_Health/wearable_copilot_openai.py` + `Consumer_Health/Wearable_Analysis/health_copilot.py` translate wearable JSON into schema-validated action plans.
*   **Clinical Ops Automator:** `Clinical/openai_clinical_ops_automator.py` emits ICD-10/CPT suggestions, SOAP notes, and prior auth packets with local JSON validation.
*   **Lab Automation Bridge:** `Lab_Automation/openai_lab_automation_bridge.py` wraps Experiment Designer outputs in payloads.
*   **Documentation:** See [OpenAI_Health_STACK.md](OpenAI_Health_STACK.md) for workflows, CLI usage, and BioKernel integration.

### Co-Worker Stack
*   **Inbox Router:** `Clinical/anthropic_inbox_router.py` fans work items into coworkers via the Event Bus.
*   **Prior Auth Coworker:** `Clinical/Prior_Authorization/appeals_agent.py` mirrors reasoning traces.
*   **Regulatory Coworker:** `Pharma/Regulatory_Affairs/anthropic_regulatory_coworker.py` drafts CTD responses with citations.
*   **Pharmacovigilance Monitor:** `Clinical/Safety/pharmacovigilance_monitor.py` triages safety signals and emits audit-ready traces.
*   **Regulatory Drafter:** `Anthropic_Health_Stack/regulatory_drafter.py` drafts regulatory submissions with audit trails.

## 📂 Directory Structure

```text
Skills/
├── 3D_Genomics/          # Hi-C analysis, chromatin conformation
├── Agentic_AI/           # Architectures (ReAct, Plan&Solve, Orchestrators)
├── Anthropic_Health_Stack/ # Regulatory drafting, care coordination
├── Clinical/             # MedPrompt, Note Summarization, Adaptive Trials, Clinical Databases
├── Computer_Science/     # Graph Algo, Distributed Systems (Async)
├── Consumer_Health/      # Wearable analysis, symptom checking
├── Data_Visualization/   # Heatmaps, volcano plots, circos, genome tracks
├── Drug_Discovery/       # ChemCrow, Self-Driving Labs, Chemoinformatics
├── Epigenomics/          # ChIP-seq, ATAC-seq, methylation, epitranscriptomics
├── Experimental_Design/  # Power analysis, sample size, batch design
├── External_Collections/ # Consolidated external skill libraries (see below)
├── Foundation_Models/    # AlphaFold3 Wrapper, BiomedGPT
├── Gene_Therapy/         # AAV vector design
├── Genomics/             # Single Cell, Spatial, CRISPR, Variant Calling, Assembly
├── Hematology/           # AML, MPN, Flow Cytometry, Coagulation
├── Imaging_Analysis/     # Imaging mass cytometry, spatial analysis
├── Immunology_Vaccines/  # CAR-T, TCR/BCR analysis, epitope prediction
├── Lab_Automation/       # Opentrons, PyLabRobot, self-driving labs
├── LLM_Research/         # RAG, Fine-Tuning, Prompt Engineering
├── Longevity_Aging/      # Biological age clocks, senescence
├── Machine_Learning/     # Biomarker discovery, survival analysis, model interpretation
├── Mathematics/          # Bayesian Opt, Linear Algebra, Probability
├── MCP_Servers/          # BioMCP Implementation
├── Metabolomics/         # Peak detection, annotation, pathway mapping
├── Microbiome/           # 16S/ITS amplicon, diversity, differential abundance
├── Multi_Omics/          # Cross-modality integration, factor analysis
├── NGS_QC/               # Read QC, alignment, trimming
├── Oncology/             # ctDNA, Liquid Biopsy, Radiomics, Pathomics
├── Pathology_AI/         # Digital pathology, computational pathology
├── Pharma/               # Drug interactions, regulatory affairs
├── Population_Genetics/  # GWAS, phylogenetics, epidemiological genomics
├── Precision_Medicine/   # Pharmacogenomics, polygenic risk scores
├── Protein_Science/      # AlphaFold3, ESM3, protein design
├── Proteomics/           # Mass spec, deep visual proteomics
├── Quantum_Biotech/      # Quantum docking, molecular simulation
├── Radiology_AI/         # Chest X-ray, CT analysis, report generation
├── Research_Tools/       # Literature mining, database access, reporting, pathway analysis
├── RNA_Therapeutics/     # ASO, siRNA, mRNA vaccine, LNP design
├── Self_Driving_Labs/    # Autonomous lab control, Bayesian experiment design
├── Sequence_Analysis/    # Alignment, SAM/BAM, FASTA/FASTQ, primer design
├── Software_Engineering/ # Web Dev (React/Next.js), Data Science (Pandas), Core Python
├── Structural_Biology/   # CryoEM, protein structure prediction
├── Synthetic_Biology/    # Genetic circuits, metabolic pathways
├── Systems_Biology/      # Flux balance analysis, metabolic reconstruction
├── Transcriptomics/      # Differential expression, RNA quantification, splicing
├── Variant_Interpretation/ # ACMG classification, clinical interpretation
└── Workflow_Management/  # Snakemake, Nextflow, CWL pipelines

```

### External_Collections
Consolidated third-party skill libraries and frameworks:
- **Auto-GPT** / **Auto-GPT-Plugins**: Autonomous GPT frameworks
- **awesome-llm-skills**: Community-curated LLM skills
- **semantic-kernel**: Microsoft's AI orchestration SDK
- **langchain**: LangChain framework components
- **google-gemini-cookbook**: Google Gemini examples
- **mcp-servers-reference**: MCP server implementations
- **Antigravity_Skills**: Universal SKILL.md agents
- **Awesome-Biomedical-LLM-Agents**: Biomedical agent resources
- **life-sciences_Claudeai-main**: Life sciences Claude integrations

### bioSkills Integration (New Feb 2026!)
Comprehensive bioinformatics skills from the bioSkills collection (388 skills):
- **Sequence_Analysis/**: FASTA/FASTQ handling, alignment, primer design
- **NGS_QC/**: FastQC, read trimming, short-read alignment
- **Transcriptomics/**: DESeq2, edgeR, alternative splicing, Ribo-seq
- **Epigenomics/**: ChIP-seq, ATAC-seq, methylation, CLIP-seq
- **3D_Genomics/**: Hi-C contact matrices, TADs, compartments
- **Population_Genetics/**: GWAS, phylogenetics, pathogen genomics
- **Data_Visualization/**: Publication-quality plots
- **Machine_Learning/**: Biomarker discovery, survival analysis
- **Multi_Omics/**: MOFA, SNF, cross-modality integration
- **Workflow_Management/**: Snakemake, Nextflow, 35+ ready pipelines
- **Systems_Biology/**: Flux balance analysis, metabolic modeling

## 🛠️ Usage Examples

**1. Run the Multi-Agent Orchestrator:**
```bash
python3 Agentic_AI/Multi_Agent_Systems/orchestrator.py
# Tutorial: Agentic_AI/Multi_Agent_Systems/TUTORIAL_BUILDING_A_SWARM.md
```

**2. Transpile a Universal Skill:**
```bash
python3 ../platform/optimizer/usdl_transpiler.py --file ../platform/optimizer/my_skill.json
# Tutorial: ../platform/optimizer/TUTORIAL_USDL_TRANSPILER.md
```

**3. Run a Self-Driving Lab Simulation:**
```bash
python3 Mathematics/Probability_Statistics/bayesian_optimization.py
# Tutorial: Mathematics/Probability_Statistics/TUTORIAL_SELF_DRIVING_LAB_OPTIMIZER.md
```

**4. Annotate Single-Cell Data:**
```bash
python3 Skills/Genomics/Single_Cell/Cell_Type_Annotation/RNA/universal_annotator.py
# Tutorial: Skills/Genomics/Single_Cell/TUTORIAL_CELL_ANNOTATION.md
```


## 📈 Roadmap (2026)
*   [x] **Phase 1:** Core Architectures (Orchestrator, Async Runtime) - *Completed Jan 2026*
*   [x] **Phase 2:** Math Foundations (Bayesian Opt, Graph Theory) - *Completed Jan 2026*
*   [x] **Phase 3:** Single Cell & Clinical Simulators - *Completed Jan 2026*
*   [x] **Phase 4:** Initial MCP Server Integration - *Completed Jan 2026*
*   [x] **Phase 5:** Dual Health Stacks (OpenAI/Anthropic) & USDL - *Completed Jan 2026*
*   [ ] **Phase 6:** Deployment to FHIR Servers & Real Lab Integration

---
*Maintained by the Artificial Intelligence Group.*


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->