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



<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->

---
name: 'precision-grounded-variant-summarization'
description: 'Summarize genetic variants with LLM assistance grounded in evidence databases, provenance, conflicts, and hallucination controls.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Precision Grounded Variant Summarization

## Overview

This skill helps produce trustworthy narrative summaries for genetic variants by grounding each claim in evidence-based databases and source documents. It is intended for clinical genomics, translational research, and variant curation contexts where free-form generation is not sufficient because provenance, uncertainty, and conflicting evidence must be visible.

## When to Use This Skill

- A user asks for a clinical or research summary of a genetic variant, gene, disease association, or variant-drug relationship.
- The task requires database-grounded evidence rather than a general explanation from model memory.
- The output must distinguish established evidence from uncertain, conflicting, missing, or low-confidence information.
- The user needs provenance for claims, including database names, access dates when available, identifiers, and source links.
- The user asks for a concise variant interpretation narrative but does not require a full ACMG/AMP classification workflow.
- The task involves hallucination control, evidence reconciliation, or audit-ready summarization for clinical genomics communication.

## Core Capabilities

1. **Variant normalization and scope definition** - Convert the user request into precise identifiers when possible, including gene symbol, transcript, HGVS expression, rsID, ClinVar Variation ID, genomic coordinates, disease context, inheritance model, and specimen or cancer context.
2. **Evidence source selection** - Prefer authoritative databases and primary sources appropriate to the question, such as ClinVar, ClinGen, OMIM, dbSNP, gnomAD, COSMIC, CIViC, OncoKB, PharmGKB, PubMed, and local files supplied by the user.
3. **Claim-grounded extraction** - Extract only claims supported by retrieved evidence, keeping the source, identifier, date, assertion type, condition, submitter or expert panel when relevant, and the exact clinical context attached to each claim.
4. **Conflict and uncertainty handling** - Identify discordant classifications, mismatched disease terms, outdated assertions, population-frequency concerns, limited evidence, and differences between germline, somatic, pharmacogenomic, and functional contexts.
5. **Narrative summarization with provenance** - Generate a concise summary that cites the evidence basis for pathogenicity, benignity, therapeutic relevance, phenotype association, functional effect, or lack of evidence.
6. **Hallucination controls** - Do not fill missing fields from memory, do not infer clinical significance from gene function alone, and label absent evidence as not found in the consulted sources rather than as evidence of absence.
7. **Audit-ready output review** - Check that every clinically meaningful statement maps to a cited source and that limitations, database disagreements, and required human review are explicit.
8. **Precision Grounding architecture** - Select authoritative databases according to the variant question and clinical context; normalize variant identifiers before cross-source joins; map retrieved records into a common structure without erasing source-specific context; and attach evidence-level provenance, including the database, record identifier, assertion or classification, condition, version or access date, and citation, to each generated claim.
9. **Precision Grounding evaluation safeguards** - Preserve and explicitly report conflicting classifications; check for stale records using available update, version, or access dates; audit citation completeness for clinically meaningful claims; test generated summaries for unsupported or hallucinated statements; and require domain-expert review of generated variant summaries and unresolved conflicts.
10. **Version-aware evidence reconciliation and abstention** - Normalize variant identity before joining records; join evidence only across explicit genome build, transcript, allele, condition, and database version or access-date context; preserve conflicting classifications instead of collapsing them; attach claim-level provenance and traceable citations; flag evidence freshness limitations; and abstain from asserting a classification or relationship when support is insufficient, ambiguous, or not traceable.

## Inputs / Outputs

**Inputs**

- Variant details: HGVS, rsID, ClinVar Variation ID, genomic coordinate, VCF row, protein change, or a free-text variant mention.
- Context: gene, disease or phenotype, inheritance pattern, germline or somatic setting, cancer type, drug context, population, and transcript/build when known.
- Evidence sources: database links, PubMed IDs, local reports, VCF annotations, curated spreadsheets, or permission to fetch public database pages.
- Output preferences: concise paragraph, structured table, patient-facing language, curator note, or clinician-facing summary.

**Outputs**

- Normalized variant identity and any unresolved ambiguity.
- Evidence table with source, identifier, assertion, condition, date or version when available, and relevance to the requested context.
- Grounded narrative summary that separates variant identity, clinical significance, disease or phenotype association, population frequency, functional evidence, therapeutic relevance, and limitations.
- Conflict report describing disagreements between sources or submitters and how they affect confidence.
- Provenance list with URLs, PubMed IDs, database records, and local file references used.
- Final caution that the summary supports expert review and is not a standalone clinical diagnosis or treatment recommendation.

## References

- Source finding: PubMed 41950627, "Precision Grounding: augmenting large language models with evidence-based databases for trustworthy genetic variant summarization." https://pubmed.ncbi.nlm.nih.gov/41950627/
- ACMG/AMP sequence variant interpretation standards and guidelines. https://pubmed.ncbi.nlm.nih.gov/25741868/
- ClinVar as a public archive for variant-phenotype relationships. https://pubmed.ncbi.nlm.nih.gov/24234437/
- CIViC knowledgebase for clinical interpretation of cancer variants. https://pubmed.ncbi.nlm.nih.gov/28138153/
