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
name: 'clinical-cost-effectiveness-modeling'
description: 'Build auditable cost-effectiveness models for clinical interventions, including AI-derived outcomes, with QALYs, ICERs, uncertainty analysis, validation, and human review.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Clinical Cost-Effectiveness Modeling

## Overview

Use this skill to translate clinical outcomes, resource use, costs, and health-state utilities into an auditable economic evaluation of competing interventions. It supports cohort and patient-level models, including analyses in which outcomes are extracted or predicted by AI, while preserving source provenance, uncertainty, validation, and human decision oversight.

## When to Use This Skill

- Compare a clinical intervention with usual care, no intervention, or another active strategy using costs and health outcomes.
- Build or review a decision tree, Markov/state-transition model, semi-Markov model, partitioned-survival model, discrete-event simulation, or patient-level microsimulation.
- Estimate quality-adjusted life-years (QALYs), incremental costs, incremental effects, incremental cost-effectiveness ratios (ICERs), or net monetary benefit (NMB).
- Extrapolate trial, registry, electronic health record, or observational evidence beyond observed follow-up.
- Incorporate AI-extracted or model-predicted clinical outcomes into a health-economic model without obscuring their provenance or error.
- Run deterministic, scenario, threshold, or probabilistic sensitivity analyses and communicate decision uncertainty.
- Calibrate a model to external targets, validate it against independent evidence, or audit an existing cost-effectiveness analysis.
- Prepare a reproducible technical report for health technology assessment, payer, provider, research, or policy review.

Do not use the model as a substitute for clinical judgment, budget-impact analysis, affordability assessment, or jurisdiction-specific health technology assessment requirements. Obtain human review before results inform care, reimbursement, procurement, or policy.

## Core Capabilities

1. **Frame the decision problem.** Define the target population, intervention, comparator, perspective, setting, time horizon, cycle length, outcome measure, currency and price year, discounting convention, willingness-to-pay context, and decision owner. Record the source and rationale for every choice.

2. **Choose the model structure.** Use the simplest structure that represents the clinical pathway and relevant timing. Apply a decision tree to short, non-recurring pathways; a cohort state-transition model to recurring population-level risks; and patient-level simulation when history, heterogeneity, competing risks, or event timing materially affects results. State whether transitions are memoryless, tunnel-state, semi-Markov, or history-dependent.

3. **Create an assumption and evidence register.** Assign each parameter and structural choice a stable identifier. Record definition, value, unit, subgroup, source, extraction method, transformation, uncertainty, correlation, date accessed, and reviewer status. Distinguish observed data, published estimates, expert elicitation, calibrated values, and unverified assumptions. Never fill evidence gaps with invented values.

4. **Specify clinical events and transitions.** Make states mutually exclusive and collectively sufficient for the modeled question. Document allowed transitions, event ordering, competing risks, recurrence, treatment switching, discontinuation, adverse events, mortality, and absorbing states. Convert rates, risks, odds, or hazards only with an explicit formula and compatible time interval; do not assume proportional hazards or constant hazards without justification.

5. **Model costs and resource use.** Include cost categories consistent with the declared perspective, such as intervention, hospitalization, monitoring, medicines, adverse events, downstream care, patient expenses, informal care, and productivity. Keep quantities separate from unit costs, use a common price year and currency, document inflation or conversion methods, and avoid double counting bundled resources.

6. **Model utilities and QALYs.** Map each health state or event to an appropriate preference-based utility, identify the instrument and population when known, and keep temporary disutilities separate from persistent state utilities. Accumulate survival-weighted utility over time, apply a justified within-cycle correction, and test any extrapolation or mapping assumption. Flag values outside plausible utility bounds or combinations that double count quality-of-life loss.

7. **Apply time and discounting consistently.** Align time zero, cycle boundaries, event timing, and the cost and outcome discounting convention. For annual discount rate `r` and time `t` in years, use `PV = value / (1 + r)^t` unless the jurisdiction requires another convention. Report undiscounted totals when they materially aid interpretation.

8. **Calculate incremental outcomes correctly.** For intervention `A` versus comparator `B`, calculate `Delta C = C_A - C_B`, `Delta E = E_A - E_B`, and `ICER = Delta C / Delta E` only when the ratio is interpretable. Identify simple and extended dominance before ranking strategies. Also report `NMB(lambda) = lambda * Delta E - Delta C`, because NMB remains interpretable when effects are near zero or the intervention is dominated.

9. **Run deterministic and structural sensitivity analyses.** Vary influential parameters over evidence-based ranges, record the range source, and report one-way, multi-way, threshold, and scenario analyses as appropriate. Test structural alternatives for model type, horizon, perspective, extrapolation, treatment waning, utility source, costing approach, adherence, and AI-derived outcome classification rather than treating them as ordinary parameter noise.

10. **Run probabilistic sensitivity analysis (PSA).** Assign distributions that respect parameter support and evidence: beta for probabilities or bounded utilities when suitable, gamma or lognormal for positive costs, Dirichlet for multinomial transitions, and multivariate draws when correlations matter. Preserve parameter dependence, use reproducible random seeds, check simulation stability, and report incremental cost-effect pairs, cost-effectiveness acceptability results, and expected NMB across relevant thresholds. Do not claim that PSA resolves structural or source uncertainty.

11. **Calibrate transparently.** Predefine calibration targets, target uncertainty, fit statistics, parameter bounds, algorithms, stopping rules, and acceptance criteria. Separate calibration targets from validation targets, retain accepted parameter sets when uncertainty propagation requires them, and show how calibration affects decision results. Reject non-identifiable or clinically implausible solutions even when numerical fit appears adequate.

12. **Validate at multiple levels.** Check formulas, units, probability sums, conservation of the cohort, impossible transitions, extreme inputs, and independently calculated test cases. Assess face validity with clinical and health-economic reviewers, internal validity against source data and expected limiting behavior, cross-validity against comparable models, and external validity against unused populations, periods, or datasets. Explain transportability limits rather than labeling a model universally valid.

13. **Govern AI-derived inputs.** Preserve the source text or record identifier, extraction prompt or model version, preprocessing, classification rule, confidence or adjudication status, and link from each derived outcome to its model parameter. Evaluate extraction or prediction error separately from economic-model uncertainty, propagate plausible misclassification scenarios, protect patient data, and require human adjudication for ambiguous or decision-critical records.

14. **Communicate uncertainty and limitations.** Present absolute and incremental discounted costs and QALYs, ICER or dominance status, NMB, key drivers, subgroup findings, sensitivity results, and model limitations. Separate statistical, parameter, structural, methodological, heterogeneity, and generalizability uncertainty. Use conditional language and avoid converting a threshold comparison into a clinical recommendation.

15. **Maintain an audit trail and human review gate.** Version the model specification, input table, code or formulas, random seed, run configuration, outputs, and change log. Require named clinical, economic, and data review as appropriate before decision use. Record reviewer comments, unresolved issues, conflicts of interest, and approval status; stop and label the analysis exploratory if required evidence or review is missing.

For a rapid execution, proceed in this order: validate the decision question and minimum inputs; build the base case; run deterministic checks; add PSA only when uncertainty inputs are defensible; perform calibration or external validation only when independent targets exist; then produce the audit package. If the available inputs cannot support a valid result within the time box, return a model plan, evidence gaps, and a clearly labeled non-decision-grade status instead of fabricating completion.

## Inputs / Outputs

### Inputs

- Decision specification: population, interventions and comparators, perspective, setting, time horizon, cycle length, outcome, price year, currency, discount rates, and threshold context.
- Clinical evidence: baseline risks, treatment effects, event timing, adverse events, discontinuation, mortality, recurrence, and subgroup modifiers, with uncertainty and provenance.
- Economic evidence: resource quantities, unit costs, utility weights, disutilities, productivity or caregiver effects when relevant, and inflation or currency-conversion sources.
- Model constraints: candidate structures, structural assumptions, extrapolation functions, transition rules, half-cycle or other within-cycle correction, and implementation environment.
- Validation evidence: internal summaries, calibration targets, independent external datasets or published targets, expert review criteria, and jurisdictional reporting requirements.
- AI-derived data when used: source-record identifiers, extraction or prediction method and version, evaluation results, missingness, confidence, adjudications, privacy constraints, and error scenarios.

Minimum viable inputs are a defined population, at least two alternatives, a perspective, a time horizon, compatible cost and effect estimates, and traceable sources. If any are absent, request them or produce only a provisional analysis plan.

### Outputs

- A decision-problem statement and model schematic or transition table.
- A machine-readable input table and human-readable evidence/assumption register with units, distributions, correlations, and provenance.
- A reproducible base-case model with discounted and undiscounted totals where appropriate.
- Incremental costs, incremental QALYs or other effects, dominance assessment, ICER when interpretable, and NMB across stated thresholds.
- Deterministic and structural sensitivity results, plus PSA outputs when valid uncertainty distributions are available.
- Calibration diagnostics and external-validation results when suitable targets exist.
- Tables or figures showing cost and QALY composition, influential assumptions, incremental cost-effectiveness, and decision uncertainty.
- A limitations and uncertainty narrative that distinguishes evidence gaps from model findings.
- An audit package containing model version, formulas or code, parameter provenance, run settings, random seed, checks, reviewer status, and unresolved issues.
- A prominent status label: `decision-grade after human review`, `exploratory`, or `insufficient evidence`, with the rationale for that label.

## References

- Kamalmaz A, Villarroya C, Mayer MA, Leis A, Garcia-Alzorriz E, et al. “Cost-effectiveness analysis of resective epilepsy surgery in drug-resistant patients: an artificial intelligence data modeling.” *Epilepsy & Behavior*. 2026;180:111027. [PubMed PMID 41946095](https://pubmed.ncbi.nlm.nih.gov/41946095/)
