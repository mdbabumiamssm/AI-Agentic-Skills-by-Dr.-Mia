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

---
name: precision-oncology-agent
description: Fuse genomic variants, pathology findings, and clinical context to draft evidence-linked therapy options for tumor board review.
measurable_outcome: Deliver a ranked therapy list with OncoKB and NCCN citations plus a data-gap checklist for every case within 10 minutes of receiving inputs.
allowed-tools:
  - read_file
  - run_shell_command
---

## At-a-Glance
- **description (10-20 chars):** Tumor board copilot
- **keywords:** oncology, genomics, OncoKB, therapy-ranking, evidence
- **measurable_outcome:** Deliver a ranked therapy list with OncoKB/NCCN citations plus data-gap checklist for every case within 10 minutes of receiving inputs.

## Inputs
- `vcf_path` (hg38 preferred) plus optional CNV/fusion summaries.
- `pathology_report` text for histology/grade/biomarkers.
- `clinical_context` dict capturing tumor type, stage, prior lines, ECOG.

## Outputs
1. Ranked treatment options (approved, off-label, clinical trials) with evidence strength + contraindications.
2. Variant interpretation table (pathogenicity, tier, therapy linkage).
3. Biomarker summary (TMB, MSI, PD-L1 if provided) and missing-test checklist.

## Core Capabilities
- Incorporate ASCO GI cancer AI management patterns by mapping each case across diagnosis, imaging/endoscopy interpretation, treatment selection, longitudinal monitoring, and real-world data inputs such as EHR and protocol-derived trial-matching context; keep every AI-generated synthesis under explicit human tumor-board governance.
- For GI cancer cases, use AI-assisted management as a specialty example by synthesizing evidence-linked use cases across diagnosis, staging, treatment selection/sequencing, radiology and pathology support, surveillance planning, trial matching/eligibility, and toxicity monitoring; route outputs through explicit clinician oversight in tumor-board workflows and state that recommendations are evidence-organizing support that must be reconciled with current guideline-concordant oncology recommendations, clinician judgment, patient preferences, and local trial availability.
- For GI cancer workflows, integrate multimodal inputs from endoscopy findings, radiology staging, pathology/histology and biomarkers, genomic results, and clinical context to draft guideline-grounded treatment options and data gaps; keep AI limited to synthesis and option drafting, with final treatment selection, sequencing, trial eligibility, and exceptions reserved for human clinician/tumor-board review.
- For GI cancer AI management summaries, explicitly cover endoscopy, radiology, pathology, molecular profiling, treatment response assessment, and trial matching; grade each evidence linkage and require tumor-board review before any clinical use rather than generating free-form treatment recommendations.

## Workflow
1. **Ingest & normalize:** Harmonize gene symbols, genome build, and variant effects.
2. **Annotate:** Query OncoKB/NCCN + internal knowledge for actionability tiers.
3. **Contextualize:** Blend pathology + prior therapy info to filter contraindicated options.
4. **Recommend:** Present therapies ordered by evidence + patient fit; cite sources.
5. **Gaps:** Highlight assays or confirmations still required before treatment.

## Guardrails
- No autonomous treatment decisions—flag outputs as advisory.
- Cite evidence rigorously (guideline version, publication).
- Highlight resistance mechanisms and prior exposure conflicts.

## References
- See `README.md` for detailed workflow plus cited Nature Cancer study.
- Harnessing Artificial Intelligence for the Management of Patients With GI Cancers. PubMed PMID: 42044465. https://pubmed.ncbi.nlm.nih.gov/42044465/


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->
