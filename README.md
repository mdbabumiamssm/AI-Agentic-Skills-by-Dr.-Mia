# Universal Biomedical Skills & Agents (2026 Edition)

![Status](https://img.shields.io/badge/Status-Active-green)
![Agents](https://img.shields.io/badge/Agents-Orchestrated-blue)
![Domain](https://img.shields.io/badge/Domain-Biotech%20%7C%20Clinical%20%7C%20Genomics-purple)
![Tech](https://img.shields.io/badge/Tech-MCP%20%7C%20DeepSeek%20%7C%20Gemini-orange)

## 🚀 Overview

This repository acts as a **Biomedical Operating System**, hosting a comprehensive library of **skills, agents, and mathematical foundations** for modern (2026) Artificial Intelligence. 

Unlike standard chatbot repos, this project focuses on **Agentic Workflows**—where autonomous systems plan, execute, use tools, and correct themselves to solve complex scientific problems. It is designed to support high-impact research, clinical decision support, and automated lab operations.

## 👤 Author & Maintainer

**MD BABU MIA, PhD**  
*Professor of Hematology & Medical Oncology | Mount Sinai*  
Specializing in MPN Research, Single-Cell Multi-Omics, and Computational Biology.

---

## 🌟 Key Capabilities (New for 2026)

### 🧬 Genomics & Single Cell
*   **Universal Annotator:** `Skills/Genomics/Single_Cell/Cell_Type_Annotation/RNA/universal_annotator.py` wraps Marker-based, Deep Learning (CellTypist), and LLM-based annotation strategies.
*   **CellAgent:** `Skills/Genomics/Single_Cell/CellAgent/SKILL.md` - Autonomous scRNA-seq analysis agent (QC -> Annotation -> DE).
*   **Variant Interpretation:** `Skills/Genomics/Variant_Interpretation/SKILL.md` - Automated ACMG classification for genetic variants.

### 🧠 Agentic AI (The Brain)
*   **Swarm Orchestrator:** `Skills/Agentic_AI/Multi_Agent_Systems/orchestrator.py` - Coordinates specialized agents (Researcher, Reviewer, Safety) to solve complex queries.
*   **Self-Correction:** Implements Reflexion patterns for iterative improvement of outputs.

### 🏥 Clinical & Operations
*   **Prior Auth Coworker:** `Skills/Clinical/Prior_Authorization/appeals_agent.py` - Drafts insurance appeals with policy-compliant reasoning traces.
*   **Clinical NLP:** `Skills/Clinical/Clinical_NLP/entity_extractor.py` - Extracts entities (Problems, Meds) from unstructured notes.
*   **Trial Matching:** Matches patient profiles to clinical trials using intelligent criteria mapping.

### 💊 Pharma & Regulatory (New!)
*   **Regulatory Drafter:** `Skills/Pharma/Regulatory_Affairs/SKILL.md` - Automates drafting of FDA CTD sections with citation management.
*   **Molecule Evolution:** `Skills/Drug_Discovery/Molecule_Design/SKILL.md` - Iteratively designs de novo drugs using LLM-driven genetic algorithms.
*   **Pharmacovigilance:** Monitors safety signals and emits audit-ready traces.

## 📂 Directory Structure

The repository is organized into domain-specific modules:

```text
Skills/
├── Agentic_AI/           # Orchestrators, Swarms, Planning Agents
├── Clinical/             # NLP, Prior Auth, Clinical Trials, FHIR
├── Drug_Discovery/       # Molecule Design, ChemCrow, Protein Structure
├── Genomics/             # Single Cell, CRISPR, Variant Interpretation
├── Pharma/               # Regulatory Affairs, Pharmacovigilance
├── Mathematics/          # Bayesian Opt, Linear Algebra
├── Software_Engineering/ # Best Practices (React, Python, Pandas)
└── External_Collections/ # Third-party skill libraries
```

## 🛠️ Usage Examples

**1. Run the Multi-Agent Orchestrator:**
```bash
python3 Skills/Agentic_AI/Multi_Agent_Systems/orchestrator.py --mission "Investigate CRISPR usage on human embryos."
```

**2. Draft a Regulatory Document:**
```bash
# Uses the new Regulatory Drafter skill
python3 Skills/Anthropic_Health_Stack/regulatory_drafter.py --input "./data/tox" --section "2.4"
```

**3. Design a Molecule:**
```bash
# Evolve a binder for a target protein
python3 Skills/Drug_Discovery/Molecule_Design/evolution_agent.py
```

## 📜 Standardized Skill Format

As of Jan 2026, all skills are being migrated to a standardized `SKILL.md` format featuring:
*   **Keywords:** For easier indexing and discovery.
*   **Measurable Outcomes:** SMART goals (e.g., "Achieves >85% accuracy").
*   **Compatibility:** explicit system and library requirements.

### 🔁 Skill Maintenance Checklist (2026+)

To keep the catalog consistent, follow this workflow whenever updating or adding skills:

1. **Scan `Skills/`:** Locate every existing `SKILL.md` (use `rg --files -g 'SKILL.md' Skills`).
2. **Enhance Each Entry:** Ensure the frontmatter or first section includes:
   - `description` (concise 10–20 characters)
   - `keywords` (3–5 core terms)
   - `measurable_outcome` (SMART, time-bound)
3. **Add New Skills:** When you find coverage gaps, scaffold a fresh skill (use `scripts/init_skill.py` if applicable) so every major workflow has a `SKILL.md`.
4. **Stage Changes:** `git add Skills/` – stage only the skill directories you touched.
5. **Commit:** Use a conventional commit message such as `feat: Enhance and add SKILL.md` (required for catalog-wide refreshes).
6. **Push:** `git push origin main` to publish the refresh so GitHub reflects the latest README/SKILL updates immediately.

Following these steps keeps the metadata discoverable and ensures reviewers can see the delta directly on GitHub (including this README).

## 📄 License

Proprietary - For personal use by MD BABU MIA.
MIT License applies to open-source components where noted.

---
*Updated: January 27, 2026*
