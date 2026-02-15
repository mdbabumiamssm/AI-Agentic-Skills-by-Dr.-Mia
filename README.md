# Universal Biomedical Skills & Agents (Biomedical OS - 2026)

![Status](https://img.shields.io/badge/Status-Active-green)
![Architecture](https://img.shields.io/badge/Architecture-Biomedical%20OS-blueviolet)
![Domain](https://img.shields.io/badge/Domain-Biotech%20%7C%20Clinical%20%7C%20Genomics-purple)
![Tech](https://img.shields.io/badge/Tech-MCP%20%7C%20DeepSeek%20%7C%20Gemini-orange)

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

## 🚀 The "Best of Best" Biomedical Skills Collection

**This repository represents the gold standard for Agentic AI in Biomedicine.**

We have painstakingly cataloged, standardized, and optimized over **600+ specialized skills** covering the entire spectrum of modern biomedical science. From single-cell genomics to clinical decision support, every tool in this library has been upgraded to a unified **Agentic Standard**.

This is not just a code dump; it is a **Biomedical Operating System (BioOS)** designed to orchestrate autonomous research.

---

## 👤 Author & Maintainer

**MD BABU MIA, PhD**  
*Assistant Professor of Hematology & Medical Oncology, Machine Learning -AI | Mount Sinai*  
Mount Sinai Tisch Cancer Institute
Icahn School of Medicine at Mount Sinai
Mount Sinai Hospital
One Gustave L. Levy Place
New York, NY 10029
Desk phone:(212) 241-2764 (x42764)
Mobile phone:(332) 256-3038
Email: md.babu.mia@mssm.edu
Specializing in Hemato-Oncology,and Machine Learning-LLM-AI.

---

## 🌟 Major Updates (February 2026)

We have successfully migrated the entire codebase to the **Universal Skill Definition Language (USDL)** standard.

### 📜 Standardized Skill Format (SKILL.md)
Every one of the 600+ skills now adheres to a rigorous metadata schema, ensuring they are **discoverable, measurable, and safe** for autonomous execution.

*   **Discoverable:** Indexed by keywords and domain.
*   **Measurable:** Every skill defines a "Measurable Outcome" (e.g., "Analyze 10k cells in <5 mins").
*   **Safe:** "Allowed Tools" restrict agent capabilities to prevent unintended side effects.

### 📚 The Catalog
We have implemented a dynamic cataloging system. The entire library is indexed in `skills_catalog.json`.

**To generate the latest catalog:**
```bash
python3 platform/skills_catalog.py
```

### 🧬 Key Domain Highlights

#### 1. Genomics & Bioinformatics
*   **BioMaster:** Orchestrator for RNA-seq, ChIP-seq, and Hi-C pipelines.
*   **CellAgent:** Autonomous single-cell annotation and quality control.
*   **STAgent:** Spatial transcriptomics analysis for Visium/Xenium.
*   **Genome Engineer:** CRISPR gRNA design and off-target prediction.

#### 2. Clinical & Operations
*   **ChatEHR:** Clinical assistant for patient records (FHIR integrated).
*   **TrialGPT:** Intelligent patient-to-trial matching and ranking.
*   **RadGPT:** Radiology report summarizer and explainer.
*   **Precision Oncology:** Multimodal treatment planning (H&E + Genomics).

#### 3. Drug Discovery
*   **MAGE:** Generative antibody design using protein language models.
*   **CheMatAgent:** Computational chemistry for molecule design.
*   **Biomni:** General-purpose biomedical research agent (150+ tools).

#### 4. Software Engineering (New!)
*   **Codebase Investigator:** Autonomous repository mapping and analysis.
*   **Technical Writing:** Automated generation of high-quality documentation.
*   **Data Visualization:** Publication-ready scientific plotting agents.

## 📂 Directory Structure

The repository is organized into high-level domains:

```text
Skills/
├── Agentic_AI/           # Orchestrators, Swarms, Planning Agents
├── Clinical/             # EHR, Radiology, Oncology, Trials
├── Drug_Discovery/       # Antibody Design, Small Molecules, Chemistry
├── Genomics/             # Single Cell, Spatial, CRISPR, Variant Interpretation
├── Hematology/           # Bone Marrow, Flow Cytometry, MPN Analysis
├── Immunology_Vaccines/  # CAR-T, TCR, Neoantigen Prediction
├── MCP_Servers/          # BioMCP and other protocol servers
├── Research_Tools/       # Biomni, Literature Mining, Knowledge Graphs
├── Software_Engineering/ # Best Practices, Codebase Analysis, Tech Writing
└── Writing_and_Productivity/ # Technical documentation and reporting
```

## 🛠️ Usage Examples

**1. Match a Patient to a Clinical Trial (TrialGPT):**
```bash
python3 Skills/Clinical/Trial_Matching/TrialGPT/run_matching.py --patient_profile ./patient.json
```

**2. Design an Antibody (MAGE):**
```bash
python3 Skills/Drug_Discovery/Antibody_Design/MAGE/generate.py --antigen "spike_protein" --count 5
```

**3. Analyze Spatial Transcriptomics (STAgent):**
```bash
python3 Skills/Genomics/Spatial_Transcriptomics/STAgent/main.py --data ./visium_data.h5ad --task "cluster_domains"
```

## 📄 License

**Copyright (c) 2026 MD BABU MIA, PhD.**  
All rights reserved.

This project is licensed under the MIT License for open-source components, but the unique architectural design and agentic workflows are the intellectual property of the author. **Attribution is mandatory.**