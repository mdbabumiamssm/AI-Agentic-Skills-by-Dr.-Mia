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
name: 'oncology-consult-survival-llm'
description: 'Guide zero-shot or fine-tuned LLM workflows for predicting cancer survival from initial oncology consultation documents with leakage control and cautious reporting.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Oncology Consult Survival LLM

## Overview

This skill guides clinically cautious workflows for using general large language models to predict cancer survival from initial oncology consultation documents. It emphasizes survival endpoint definition, zero-shot versus fine-tuned evaluation, data leakage controls, calibration, uncertainty, and reporting that is appropriate for research or decision support rather than standalone clinical judgment.

Use this skill when the task involves oncology-document prognostic modeling and the user needs a defensible workflow grounded in the 2026 PubMed-indexed finding on fine-tuning versus zero-shot LLM prediction of cancer survival from initial consultation notes.

## When to Use This Skill

- Predicting survival, mortality, or time-to-event outcomes from initial oncology consultation documents.
- Comparing zero-shot prompting with fine-tuned LLM approaches for oncology prognosis tasks.
- Designing training, validation, or test splits for clinical-note survival prediction.
- Reviewing oncology prognostic modeling plans for leakage, endpoint ambiguity, label quality, or calibration gaps.
- Producing clinically cautious summaries of LLM-derived survival risk predictions.
- Auditing whether consultation-note text contains future information that would invalidate a baseline prediction task.

## Core Capabilities

1. **Survival endpoint framing**: Define the target outcome before modeling, including prediction time zero, follow-up window, event type, censoring rule, and whether the endpoint is binary survival status or time-to-event risk.

2. **Document cohort preparation**: Confirm that each input is an initial oncology consultation document and preserve metadata needed for temporal ordering, cancer type, disease stage, treatment intent, and follow-up ascertainment.

3. **Leakage control**: Exclude future information such as post-consultation treatments, later imaging, later pathology amendments, recurrence notes, death documentation, hospice enrollment after baseline, or follow-up outcomes embedded in copied text.

4. **Zero-shot workflow design**: Build explicit prompts that ask for structured prognostic outputs, cite only baseline document evidence, separate unknown from absent information, and avoid forcing a risk estimate when the note is insufficient.

5. **Fine-tuning workflow design**: Prepare de-identified training examples with stable labels, patient-level splits, development and held-out test sets, locked preprocessing, and no overlap across documents from the same patient.

6. **Evaluation and calibration**: Report discrimination, calibration, threshold behavior, uncertainty, and clinically relevant subgroup checks without inventing or overstating performance. Prefer calibration plots, Brier score, concordance-style summaries, or sensitivity and specificity when appropriate to the endpoint.

7. **Clinical caution and reporting**: State that outputs are model-derived prognostic estimates requiring clinician review and local validation. Avoid treatment recommendations unless the user separately requests guideline-grounded clinical decision support.

8. **Comparison interpretation**: When comparing zero-shot and fine-tuned models, separate performance from operational tradeoffs such as data requirements, reproducibility, site-specific adaptation, privacy burden, and robustness to documentation style.

## Inputs / Outputs

**Inputs**

- Initial oncology consultation documents or extracted consultation-note text.
- Cohort metadata such as patient identifier, consultation date, cancer diagnosis, stage, treatment intent, and source system.
- Survival labels or follow-up data with event date, censoring date, last known alive date, and endpoint definition.
- Modeling mode: zero-shot prompting, fine-tuning, or comparative evaluation.
- Constraints for de-identification, governance, protected health information handling, and permitted tools.

**Outputs**

- A structured workflow for survival prediction from baseline oncology consultation text.
- A leakage audit checklist tailored to the available documents and labels.
- Prompt templates or fine-tuning data schema suitable for the specified endpoint.
- Evaluation plan covering discrimination, calibration, uncertainty, and subgroup review.
- Clinically cautious reporting language that describes model outputs as research or decision-support signals, not definitive predictions.

## References

- Phaterpekar T, Zeng Z, Mali Y, Leung B, Ho C. "Investigating fine-tuning versus zero-shot learning for general large language models when predicting cancer survival from initial oncology consultation documents." PubMed: https://pubmed.ncbi.nlm.nih.gov/42004490/
