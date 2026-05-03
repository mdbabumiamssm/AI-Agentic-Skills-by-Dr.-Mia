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
name: 'mutual-information-fairness-auditor'
description: 'Audit and mitigate intersectional, multiclass model fairness by estimating mutual information between prediction-derived variables and sensitive attributes.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Mutual Information Fairness Auditor

## Overview

Use this skill to evaluate model fairness through a mutual-information lens, especially when sensitive attributes are multiclass, intersecting, or structured across multiple subgroup definitions. The workflow operationalizes the MIFair idea that group fairness can be assessed by measuring statistical dependence between prediction-derived variables and sensitive attributes, then reduced through regularization or governance recommendations.

This skill is intended for legal, compliance, model-risk, and AI governance reviews where conventional single-attribute or binary fairness metrics are too narrow.

## When to Use This Skill

- A user asks for a fairness audit of a classifier, scoring model, or decision-support system.
- The audit involves intersectional groups such as race-by-gender, age-by-disability, or other combined sensitive attributes.
- The model has multiclass labels, multiclass predictions, ranked outputs, or several protected attribute categories.
- Existing fairness checks are limited to binary demographic parity, equalized odds, or single-sensitive-attribute comparisons.
- The user needs bias assessment evidence for AI governance, procurement review, model documentation, or legal risk analysis.
- The user asks how to mitigate observed group dependence during training without changing the model objective entirely.

## Core Capabilities

1. **Fairness target selection** - Map the review question to a prediction-derived variable, such as predicted class, score bin, error indicator, true-positive indicator, or false-positive indicator.
2. **Sensitive attribute construction** - Represent protected attributes individually and intersectionally, preserving multiclass and multi-attribute subgroup structures when data support them.
3. **Mutual-information audit design** - Estimate dependence between the selected prediction-derived variable and sensitive attributes, treating lower mutual information as evidence of weaker statistical dependence.
4. **Fairness-notion alignment** - Relate the audit target to common group fairness goals such as independence or separation when the selected variables correspond to those concepts.
5. **Subgroup diagnostics** - Identify which attributes, intersections, or outcome classes contribute most to dependence, while flagging sparse subgroups and unstable estimates.
6. **Mitigation guidance** - Recommend in-processing regularization, inspired by prejudice-remover style methods, that penalizes mutual information between prediction-derived variables and sensitive attributes.
7. **Governance reporting** - Produce a concise audit report with assumptions, data limitations, metric definitions, subgroup coverage, findings, and recommended next actions.

## Inputs / Outputs

### Inputs

- Dataset schema with target label, prediction outputs, model scores if available, and sensitive attributes.
- Task type, especially binary or multiclass classification.
- Policy context describing which fairness notion matters: independence, separation, error-rate parity, or another group fairness target.
- Subgroup definitions, including required intersectional combinations and any legally or institutionally protected classes.
- Optional training configuration if mitigation recommendations or in-processing regularization are requested.

### Outputs

- A mutual-information fairness audit plan specifying variables, subgroup construction, estimators, and validation checks.
- A metric table covering overall dependence, per-attribute dependence, intersectional dependence, and class- or error-specific dependence when applicable.
- A subgroup reliability note identifying small samples, missing sensitive attributes, rare classes, and estimates that should not be overinterpreted.
- A mitigation recommendation describing whether to use mutual-information regularization, adjust subgroup definitions, collect more data, or perform additional human review.
- A governance-ready summary with limitations, residual fairness risk, and traceable references.

## References

- MIFair: A Mutual-Information Framework for Intersectionality and Multiclass Fairness: http://arxiv.org/abs/2604.28030v1
