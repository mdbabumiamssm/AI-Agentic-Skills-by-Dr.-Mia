# Universal Biomedical Skills Platform (v2026.3 "Enterprise Edition")

**The Professional Agentic Operating System for Life Sciences**

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Status: Enterprise Ready](https://img.shields.io/badge/Status-Enterprise_Ready-blue)]()
[![Architecture: Event Driven](https://img.shields.io/badge/Arch-Event_Driven-green)]()
[![Platform: Multi-LLM](https://img.shields.io/badge/Platform-Claude%20%7C%20ChatGPT%20%7C%20Gemini-purple.svg)]()

---

## 👨‍💻 Author & Maintainer

**MD BABU MIA, PHD**
Assistant Professor
Mount Sinai Tisch Cancer Institute
Icahn School of Medicine at Mount Sinai
Mount Sinai Hospital
One Gustave L. Levy Place
New York, NY 10029
Desk phone:(212) 241-2764 (x42764)
Mobile phone:(332) 256-3038
Email: md.babu.mia@mssm.edu

---

## 🏢 Mission

We have matured the platform from a research prototype into a **Biomedical Agentic OS**. Our mission is to standardize how AI agents—whether running on OpenAI, Anthropic, or open models—collaborate to solve complex scientific problems.

The system now features **Swarm Orchestration** (multi-agent teaming), **Dual-Stack Alignment** (OpenAI JSON + Anthropic Thinking Blocks), and **Self-Driving Lab** capabilities.

---

## 🚀 Key Capabilities (2026 Enterprise Update)

1.  **Swarm Orchestrator:** A new `SwarmOrchestrator` (`Skills/Agentic_AI/Multi_Agent_Systems`) allows specialized agents (Researcher, Reviewer, Safety Officer) to work in parallel.
2.  **Universal Adapter (Write Once, Run Anywhere):**
    *   **Google Antigravity Support:** Native support for the `SKILL.md` format.
    *   **Multi-Provider Transpiler:** Automatically converts `SKILL.md` instructions into **OpenAI** JSON Schemas, **Anthropic** XML Thinking Blocks, or **Gemini** prompts.
3.  **Dual-Stack Health:** Native support for:
    *   **OpenAI Health Stack:** JSON Schema-enforced agents for clinical triage and ops.
    *   **Anthropic Co-Worker Stack:** "Thinking Block" agents for high-stakes regulatory and safety reasoning.
4.  **Self-Driving Labs:** Integrated Bayesian Optimization (`Skills/Mathematics/Probability_Statistics`) for autonomous experiment design.
5.  **MCP Server:** A standard JSON-RPC 2.0 server to expose all skills to tools like Claude Desktop or your own agent IDE.

---

## 📂 Comprehensive Skills Overview

### 🧠 Agentic AI & Core (New)

| Skill | Description | Key Tools |
|-------|-------------|-----------|
| **[Antigravity Skills](Skills/External_Collections/Antigravity_Skills/)** | **(NEW)** Universal `SKILL.md` agents | MCP, Universal Transpiler |
| **[Swarm Orchestrator](Skills/Agentic_AI/Multi_Agent_Systems/)** | **(NEW)** Coordinate multi-agent teams | AsyncIO, Task Routing |
| **[USDL Transpiler](platform/optimizer/)** | **(NEW)** Cross-compile skills for any model | CLI, JSON/XML Gen |
| **[BioKnowledge Graph](Skills/Computer_Science/Graph_Algorithms/)** | **(NEW)** GraphRAG for drug pathways | NetworkX, Dijkstra |

### ⚕️ Clinical Skills

| Skill | Description | Key Tools |
|-------|-------------|-----------|
| **[Care Copilot](Skills/Consumer_Health/)** | **(UPDATED)** Patient triage with agentic self-correction | JSON Mode, Self-Correction |
| **[Regulatory Drafter](Skills/Anthropic_Health_Stack/)** | **(UPDATED)** Regulatory drafting with audit trails | Chain-of-Thought, XML |
| **[Prior Auth Agent](Skills/Clinical/Prior_Authorization/)** | **(NEW)** Automated medical necessity & **Appeals** | Self-Correction, MedPrompt |
| **[Recruitment Agent](Skills/Clinical/Clinical_Trials/)** | Vector-based patient-trial matching | Vector DB, Cosine Similarity |

### 🧬 Genomics Skills

| Skill | Description | Key Tools |
|-------|-------------|-----------|
| **[CRISPR Design Agent](Skills/Genomics/CRISPR_Design_Agent/)** | **(UPDATED)** Automated gRNA design | Off-target analysis |
| **[Variant Classifier](Skills/Genomics/Variant_Interpretation/)** | **(UPDATED)** ACMG classification & AI Reporting | Prompt Optimization |
| **[Single Cell Multi-Omics](Skills/Genomics/Single_Cell/)** | Integrated analysis of RNA, Protein, ATAC | scverse, Seurat |
| **[Spatial Transcriptomics](Skills/Genomics/Spatial_Transcriptomics/)** | Multimodal spatial analysis | Dynamic code generation |

### 🧪 Lab & Drug Discovery Skills

| Skill | Description | Key Tools |
|-------|-------------|-----------|
| **[Opentrons Generator](Skills/Lab_Automation/Opentrons_Agent/)** | **(NEW)** High-level intent to robot code | Opentrons API, Meta-Prompt |
| **[Self-Driving Lab](Skills/Mathematics/Probability_Statistics/)** | **(NEW)** Bayesian Optimization for experiments | Gaussian Processes, UCB |
| **[Molecule Evolution Agent](Skills/Drug_Discovery/Molecule_Design/)** | Evolutionary de novo drug design | Genetic Algorithm |
| **[Protein Structure Prediction](Skills/Drug_Discovery/Protein_Structure/)** | AI-powered structure prediction | AlphaFold 2/3 |

---

## 🤝 Dual Healthcare Alignments

| Stack | Focus | Entry Points |
| --- | --- | --- |
| **OpenAI Health Stack** | Clinical Ops, Triage, Consumer Coaching. Uses strictly enforced JSON schemas. | `Skills/OpenAI_Health_Stack/TUTORIAL_OPENAI_HEALTH.md` |
| **Anthropic Health Stack** | Regulatory Affairs, Safety Monitoring. Uses "Thinking Blocks" for auditability. | `Skills/Anthropic_Health_Stack/TUTORIAL_ANTHROPIC_WORKFLOWS.md` |

---

## 🚀 Quick Start

### 1. Run the MCP Server (Recommended)
Connect your favorite agent client (Claude Desktop, etc.) to the BioKernel:
```bash
python3 platform/biokernel/mcp_server.py
```

### 2. Run the Swarm Orchestrator
Coordinate a multi-agent research mission:
```bash
python3 Skills/Agentic_AI/Multi_Agent_Systems/orchestrator.py
```

### 3. Transpile a Skill (CLI)
Convert a generic `SKILL.md` into an OpenAI prompt:
```bash
# Requires python script update to point to SKILL.md
python3 platform/optimizer/usdl_transpiler.py --file my_skill.json --provider openai
```

### 4. Run the Dashboard
Visualize workflows in real-time:
```bash
python3 platform/dashboard.py
```

---

## License

[MIT License](LICENSE) - Free for academic and commercial use.

## Citation

If you use these skills in your research, please cite:

```bibtex
@software{universal_life_science_skills,
  author = {Mia, MD Babu},
  title = {Universal Life Science and Clinical Skills for LLM Agents},
  year = {2026},
  publisher = {GitHub},
  url = {https://github.com/ARTIFICIALINTELLIGENCEGROUP/skills}
}
```
