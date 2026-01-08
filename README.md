# LLMs Universal Life Science & Clinical Skills

MD BABU MIA, PHD,
Assistant Professor
Mount Sinai Tisch Cancer Institute
Icahn School of Medicine at Mount Sinai
Mount Sinai Hospital
One Gustave L. Levy Place
New York, NY 10029
Desk phone:(212) 241-2764 (x42764)
Mobile phone:(332) 256-3038
Email: md.babu.mia@mssm.edu
 

**The Open-Source Operating System for Biomedical AI Agents**

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Status: Production](https://img.shields.io/badge/Status-Production-brightgreen.svg)]()
[![Platform: Multi-LLM](https://img.shields.io/badge/Platform-Claude%20%7C%20ChatGPT%20%7C%20Gemini-purple.svg)]()
[![Skills: 35+](https://img.shields.io/badge/Skills-35%2B-orange.svg)]()

---

## Mission

We build **production-ready, platform-agnostic biomedical AI skills** that empower researchers, clinicians, and developers to deploy advanced AI capabilities across any LLM interface. Whether you use **Claude**, **ChatGPT**, **Gemini**, or custom open-source models, our standardized skills deliver reproducible, validated results for real-world biomedical workflows.

---

## Skills Catalog

### Agentic AI (New Core)
*The "Brain" of the system.*

| Skill | ID | Description |
|-------|-----|-------------|
| **Plan-and-Solve Agent** | `agentic.architectures.planner` | Multi-step reasoning with planning phase |
| **Supervisor Orchestrator** | `agentic.multi_agent.orchestrator` | LangGraph-style state machine |
| **ReAct Agent** | `agentic.architectures.react` | Thought-Action-Observation loops |
| **Tree of Thought** | `agentic.reasoning.tot` | Advanced problem solving |

### Computer Science & Math
*The "Engine" of the system.*

| Skill | ID | Description |
|-------|-----|-------------|
| **Knowledge Graph** | `cs.algorithms.knowledge_graph` | Biomedical pathfinding (BFS/DFS) |
| **Bayesian Optimization** | `math.prob_stats.bayes_opt` | Experimental design for Self-Driving Labs |
| **Advanced RAG** | `llm.rag.advanced` | HyDE & Contextual Reranking |

### Clinical Skills (10 skills)

| Skill | ID | Description |
|-------|-----|-------------|
| **MedPrompt Engine** | `biomedical.clinical.medprompt` | Few-shot CoT ensemble (SOTA) |
| **Hallucination Detector** | `biomedical.clinical.safety` | CHECK framework validation |
| **Clinical Note Summarization** | `biomedical.clinical.note_summarization` | SOAP format clinical notes |
| **Trial Eligibility Agent** | `biomedical.clinical.trial_eligibility` | Patient-trial matching |
| **Precision Oncology Agent** | `biomedical.clinical.precision_oncology` | Multimodal cancer recommendations |
| **Medical Imaging AI** | `biomedical.clinical.medical_imaging` | CT/MRI/X-ray analysis (MONAI) |
| **EHR/FHIR Integration** | `biomedical.clinical.ehr_fhir_integration` | Healthcare interoperability |

### Drug Discovery Skills (6 skills)

| Skill | ID | Description |
|-------|-----|-------------|
| **Real RDKit Tools** | `biomedical.drug_discovery.chem_tools` | LogP, QED, Lipinski, PAINS, SA_Score |
| **AgentD Drug Discovery** | `biomedical.drug_discovery.agentd` | Early-stage discovery |
| **Protein Structure** | `biomedical.drug_discovery.protein_structure` | AlphaFold 2/3, OpenFold |
| **Antibody Design** | `biomedical.drug_discovery.antibody_design` | MAGE multistate design |

### Genomics Skills (7 skills)

| Skill | ID | Description |
|-------|-----|-------------|
| **scGPT Agent** | `biomedical.genomics.scgpt` | Single-cell foundation model (33M cells) |
| **Single-Cell RNA QC** | `biomedical.genomics.single_cell_qc` | MAD-based adaptive filtering |
| **CRISPR Design Agent** | `biomedical.genomics.crispr_design` | Guide RNA design & off-target |
| **Variant Annotation** | `biomedical.genomics.variant_annotation` | VEP, ClinVar, ACMG classification |

---

## Repository Structure

```
skills/
├── Skills/                           # Core Source Code
│   ├── Agentic_AI/                   # [NEW] Orchestration & Reasoning
│   │   ├── Agent_Architectures/      # Plan-and-Solve, ReAct
│   │   ├── Multi_Agent_Systems/      # Supervisor, Orchestrator
│   │   └── Reasoning_Models/         # Tree of Thought
│   │
│   ├── Computer_Science/             # [NEW] Algorithms
│   │   └── Graph_Algorithms/         # Knowledge Graphs
│   │
│   ├── Mathematics/                  # [NEW] Probability & Stats
│   │   └── Probability_Statistics/   # Bayesian Optimization
│   │
│   ├── Clinical/                     # Healthcare AI
│   │   ├── Clinical_Note_Summarization/ # MedPrompt
│   │   ├── Hallucination_Detection/     # CHECK Framework
│   │   └── ...
│   │
│   ├── Drug_Discovery/               # Cheminformatics
│   │   ├── ChemCrow_Tools/           # RDKit Integration
│   │   └── ...
│   │
│   ├── Foundation_Models/            # SOTA Models
│   │   ├── scGPT_Agent/
│   │   └── ...
│   │
│   ├── LLM_Research/                 # RAG & Prompting
│   │   └── RAG_Systems/              # HyDE, Reranking
│   │
│   ├── MCP_Servers/                  # Model Context Protocol
│   └── Research_Tools/               # General tools
│
├── test_demonstration/               # Validation suite
├── platform/                         # USDL Platform Prototype
└── docs/                             # Documentation
```

---

## Why This Repository?

| Challenge | Our Solution |
|-----------|--------------|
| Biomedical AI tools are fragmented | **Universal Skill Definition Language (USDL)** compiles once, deploys everywhere |
| AI prompts lack scientific validation | Every skill follows **peer-reviewed methodologies** with citations |
| Integration is complex | **Drop-in Python modules** work with LangChain, AutoGen, Semantic Kernel |
| Results are non-reproducible | **Statistical rigor** (MAD-based filtering, etc.) ensures consistency |
| Clinical data is siloed | **FHIR/EHR integration** enables healthcare interoperability |
| Drug discovery is slow | **Knowledge graphs & AI** accelerate target identification |
| Protein structures unknown | **AlphaFold integration** provides structure predictions |

---

## Latest Updates (Phase 2.5 Expansion)

**Computer Science & Math:**
- **Knowledge Graphs**: Implemented robust BFS/DFS pathfinding for biomedical networks.
- **Bayesian Optimization**: Added mathematical engine for Self-Driving Labs (Gaussian Processes).

**Agentic AI:**
- **Plan-and-Solve**: New agent architecture that separates planning from execution for complex reasoning.
- **Orchestrator**: Production-grade supervisor pattern for multi-agent coordination.

**LLM Research:**
- **Advanced RAG**: Implemented HyDE (Hypothetical Document Embeddings) and Contextual Reranking.

**Clinical & Drug Discovery:**
- **MedPrompt**: Upgrade to full Chain-of-Thought ensemble.
- **Real RDKit**: Replaced mocks with actual cheminformatics calculations.

---

## Author & Maintainer

**MD BABU MIA**
*Artificial Intelligence Group*
*Icahn School of Medicine at Mount Sinai*
md.babu.mia@mssm.edu

---

## License

[MIT License](LICENSE) - Free for academic and commercial use.

---

## Citation

```bibtex
@software{universal_life_science_skills,
  author = {Mia, MD Babu},
  title = {Universal Life Science and Clinical Skills for LLM Agents},
  year = {2025},
  publisher = {GitHub},
  url = {https://github.com/ARTIFICIALINTELLIGENCEGROUP/skills}
}
```
