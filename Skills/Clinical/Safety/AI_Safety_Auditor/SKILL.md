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
name: ai-safety-auditor
description: Validates clinical AI outputs for safety, bias, and hallucination risks before delivery to end-users or clinicians.
keywords:
  - ai-safety
  - compliance
  - bias-detection
  - hallucination-check
  - clinical-validation
measurable_outcome: Identifies 100% of critical safety violations and potential PHI leakage in generated clinical notes.
license: MIT
metadata:
  author: Biomedical AI Team
  version: "1.0.0"
compatibility:
  - system: Python 3.10+
allowed-tools:
  - run_shell_command
  - read_file
---

# AI Safety Auditor

The **AI Safety Auditor** is a critical "human-in-the-loop" simulator and automated guardrail system. It intercepts outputs from other clinical agents to ensure they meet medical safety standards, do not contain Protected Health Information (PHI) where inappropriate, and are free from harmful hallucinations.

## When to Use This Skill

*   As a final check before any clinical agent output is shown to a user.
*   To audit historical logs of agent interactions for compliance.
*   When detecting potential bias in diagnosis or treatment recommendations.
*   To verify that citations in a generated report actually exist (hallucination check).

## Core Capabilities

1.  **PHI Scrubbing Verification**: Ensures no identifiers leaked into non-secure outputs.
2.  **Hallucination Detection**: Cross-references generated claims against trusted knowledge bases.
3.  **Bias Scanning**: Checks for demographic or socioeconomic bias in clinical reasoning.
4.  **Contraindication Check**: Verifies treatment recommendations against patient allergies/conditions.
5.  **Clinician Override Preference Auditing**: Treats clinician overrides of clinical AI recommendations as implicit preference signals by categorizing override types, separating execution capability from alignment capability, avoiding suppression bias against correct-but-difficult recommendations, and tying reward-model updates to longitudinal outcomes in value-based care.
6.  **AI-Transformed Patient-Centered Documentation Safeguards**: Applies a pre/post error taxonomy for transformed patient-centered notes, checks for contradictions against the clinician-authored source, detects omitted clinically meaningful details, triggers clinician review for unsafe or meaning-altering changes, and treats mental health documentation with heightened sensitivity to tone, attribution, and patient meaning.
7.  **On-Prem Clinical LLM Deployment Audits**: Evaluates locally hosted clinical diagnosis workflows using distilled open-source reasoning models such as DeepSeek-R1 variants by checking benchmark drift versus parent models, diagnostic task stratification, infrastructure and privacy constraints, human oversight gates, and failure modes unique to on-premises clinical diagnosis systems.
8.  **Distilled Open Reasoning Model Deployment Checklist**: Requires on-premises clinical diagnosis deployments to document diagnosis benchmark selection, model-size/performance tradeoffs, calibration and refusal behavior, privacy constraints, and comparison against curated clinical references before clinical use.
9.  **Clinical Practice Guideline Drafting Evaluation**: Audits LLM-assisted real-time guideline development workflows with a checklist for citation traceability, PICO framing, conflict resolution, evidence grading, expert adjudication, and audit logs across drafting and review steps.
10. **LLM-Assisted Guideline Draft Checklist**: Evaluates clinical practice guideline drafts generated with LLM assistance by requiring evidence traceability, recommendation grading review, conflict-of-interest review, omission checks against intended guideline scope, and real-time human editorial oversight before draft acceptance.
11. **Fine-Grained Medical Q&A Dataset Evaluation**: Reviews domain-specific medical Q&A datasets before trustworthy medical LLM benchmarking by applying item-level correctness, omission, and harm-potential ratings; stratifying results by medical domain or topic; checking rubric reliability across reviewers; and requiring dataset documentation for source provenance, coverage, and known limitations.
12. **Test-Time Knowledge Acquisition Safety Pattern**: Requires medical decision support agents to retrieve relevant evidence before reasoning, vet source quality, surface and handle contradictions, abstain when evidence is weak, and keep retrieved evidence clearly separated from model-generated inference.
13. **Medical Decision Support Retrieval Safety Gates**: Audits test-time knowledge acquisition by requiring explicit retrieval provenance, source freshness checks, contradiction handling, predefined abstention criteria, and post-retrieval clinical safety review before recommendations.
14. **Test-Time Knowledge Acquisition Decision-Support Checks**: Validates medical decision support workflows that acquire knowledge at inference time by requiring retrieval/source validation, knowledge freshness checks, contradiction escalation, evidence provenance, and fail-closed abstention or clinician review when supporting evidence is weak.
15. **Pre-Reasoning Test-Time Knowledge Acquisition Audit**: Requires clinical agents to retrieve current vetted references before medical reasoning, log source provenance, compare pre- and post-retrieval answer changes, and escalate uncertainty when retrieved evidence conflicts.
16. **Knowledge-Augmented Recommendation Release Gate**: Before allowing test-time knowledge-augmented clinical decision-support recommendations, requires source vetting, retrieval provenance, conflict detection, hallucination checks, and evaluation against no-retrieval baselines.
17. **Medical Answer Safety Adjudication**: Scores answer correctness, clinically important omissions, and risk of harm through blinded expert review; applies severity weighting, resolves reviewer disagreements through adjudication, and retains representative unsafe-response examples for audit and remediation.

## Workflow

1.  **Intercept**: Receive candidate response from a Clinical Agent.
2.  **Scan**: Run parallel safety checks (PHI, Bias, Factuality).
3.  **Verdict**: Pass, Flag for Review, or Reject.
4.  **Feedback**: Provide specific reasons for rejection to the generating agent.

## On-Prem Clinical LLM Deployment Evaluation

For distilled open reasoning models used in local medical diagnosis systems, audit deployment readiness by:

*   Comparing diagnostic performance against the parent model and approved clinical workflow to detect benchmark drift before clinical release.
*   Selecting diagnosis benchmarks that match intended clinical diagnosis tasks and stratifying evaluation so broad aggregate scores do not hide task-specific failures.
*   Comparing model-size/performance tradeoffs before selecting a distilled model for local deployment.
*   Verifying calibration, uncertainty handling, and refusal or abstention behavior for unsafe, unsupported, or out-of-scope diagnostic requests.
*   Verifying local infrastructure, access control, logging, update, privacy, and PHI handling constraints required for on-premises use.
*   Requiring comparison against curated clinical references and clinician oversight gates before clinical use of diagnostic suggestions, escalation, abstention, and high-risk or uncertain cases.
*   Recording failure modes observed during local diagnosis testing, with special attention to deployment-specific risks created by locally hosted models.

## Example Usage

**User**: "Audit this generated discharge summary for safety."

**Agent Action**:
```bash
python3 Skills/Clinical/Safety/AI_Safety_Auditor/audit_output.py --input discharge_summary.txt --checks "all"
```

## References

*   http://arxiv.org/abs/2604.28010v1
*   https://pubmed.ncbi.nlm.nih.gov/42054574/
*   https://pubmed.ncbi.nlm.nih.gov/42062641/
*   https://pubmed.ncbi.nlm.nih.gov/42042855/
*   https://pubmed.ncbi.nlm.nih.gov/42039929/
*   https://pubmed.ncbi.nlm.nih.gov/41953846/
*   https://pubmed.ncbi.nlm.nih.gov/41908501/


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->
