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
name: 'clinical-llm-patient-education'
description: 'Translate radiology reports and clinical documents into personalized patient education while preserving meaning, uncertainty, next steps, and clinician review requirements.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Clinical LLM Patient Education

## Overview

Use this skill to convert radiology reports and other clinical documents into understandable, personalized patient education without changing their diagnostic meaning. Apply a governed workflow that explains uncertainty, identifies next steps, checks for harmful simplification, and reserves consequential guidance for professional review.

## When to Use This Skill

- Explain an MRI, CT, ultrasound, X-ray, pathology, laboratory, discharge, or consultation document to a patient or caregiver.
- Rewrite clinical information at a requested reading level or in plain language.
- Personalize an explanation using supplied clinical context, preferences, language needs, or accessibility needs.
- Compare an LLM-generated explanation with a clinician-authored explanation or approved patient-education source.
- Review patient-facing clinical text for omissions, overstatement, false reassurance, alarmist language, or unsafe recommendations.
- Prepare a draft that a qualified healthcare professional will review before delivery.

## Core Capabilities

1. **Source-grounded extraction**
   - Identify the document type, body region or clinical topic, key findings, impression, recommendations, comparison studies, limitations, and explicit uncertainty.
   - Distinguish statements present in the source from contextual explanations added for patient understanding.

2. **Meaning-preserving translation**
   - Translate technical language into plain language while retaining laterality, location, severity, timing, negation, measurements, and diagnostic qualifiers.
   - Define important terms briefly and avoid substituting a more certain diagnosis than the source supports.

3. **Personalized explanation**
   - Adapt vocabulary, length, structure, reading level, and examples to the stated audience.
   - Use only supplied patient context. Do not infer diagnoses, prognosis, identity, literacy, preferences, or treatment suitability from missing information.

4. **Uncertainty communication**
   - Preserve qualifiers such as possible, probable, indeterminate, cannot exclude, stable, and unchanged.
   - Explain what is known, what remains uncertain, and what information or follow-up may resolve that uncertainty.

5. **Next-step framing**
   - Separate source-stated recommendations from general discussion prompts.
   - Present questions the patient may ask the treating professional without prescribing treatment, changing medication, or declaring that follow-up is unnecessary.

6. **Clinician comparison**
   - When a clinician-authored explanation is available, compare coverage of findings, uncertainty, urgency, next steps, tone, and patient relevance.
   - Treat disagreements as review items rather than assuming either explanation is correct without checking the source document.

7. **MRI-report education benchmarking**
   - For MRI-report patient education evaluations, compare generated explanations against expert interpretations when available, preserving uncertainty and follow-up recommendations from the report.
   - Measure omissions and overstatements as review findings, and require clinician review before any generated explanation is delivered to a patient.

8. **Harmful-simplification detection**
   - Flag omitted critical findings, lost negation, incorrect anatomy or laterality, changed severity, unsupported causal claims, false reassurance, unnecessary alarm, and invented recommendations.
   - Check whether simplification hides incidental findings, limitations, or conditional follow-up language that could affect care.

9. **Professional-review governance**
   - Require qualified professional review before patient delivery when content may affect diagnosis, prognosis, urgency, treatment, medication, procedures, surveillance, or decisions to seek or defer care.
   - Escalate immediately when the source or context indicates a potentially urgent finding, acute deterioration, or emergency symptoms. Do not replace local emergency instructions or clinical judgment.

10. **Privacy-conscious handling**
   - Use the minimum necessary patient information and omit direct identifiers from outputs unless explicitly required in an authorized workflow.
   - Do not retrieve or combine external patient data without authorization.

11. **Quality-control workflow**
    - Complete the following sequence:
      1. Parse the source and list clinically material facts.
      2. Draft a faithful plain-language explanation.
      3. Verify every material statement against the source.
      4. Check uncertainty, negation, severity, laterality, measurements, and recommendations.
      5. Compare with clinician-authored content when provided.
      6. Run the harmful-simplification and escalation checks.
      7. Label unresolved questions and required professional review.

12. **Perioperative patient-communication workflows**
    - Provide plain-language preoperative explanations and set expectations using clinician-provided information about the procedure purpose, preparation, expected course, risks, benefits, alternatives, and patient questions.
    - Support consent discussions by organizing and explaining the supplied information, but do not conduct or document informed consent, assess decision-making capacity, or replace the surgeon's consent discussion.
    - For postoperative instructions, organize only source-stated wound or device care, activity, diet, medication, follow-up, and warning-sign guidance without adding or changing clinical directions.
    - Check readability against the stated audience and target reading level, preserve uncertainty and conditional language, and adapt presentation to supplied language, sensory, cognitive, literacy, or other accessibility needs.
    - Include teach-back prompts asking the patient to explain key preparation, recovery, and escalation instructions in their own words, and require treating-surgeon approval before patient delivery or use in consent or care decisions.

## Inputs / Outputs

### Inputs

- The complete clinical document or the relevant excerpt, including its impression or conclusion when available.
- Document type and clinical context needed to interpret terminology.
- Intended audience, preferred language, target reading level, desired length, and accessibility needs.
- Clinician-authored explanation, approved education material, or institutional guidance when comparison is requested.
- Applicable review, privacy, consent, and escalation requirements.

Do not proceed as though the record is complete when pages, sections, or referenced prior studies are missing. State the limitation and request the missing material when it could change the explanation.

### Outputs

Produce a structured draft containing:

1. **Plain-language summary**: A concise account of the main findings without adding a diagnosis.
2. **What the findings mean**: Definitions and context tied directly to the source.
3. **Uncertainty and limitations**: Explicit qualifiers, unresolved questions, and study limitations.
4. **Next steps in the source**: Recommendations stated by the report or clinical document.
5. **Questions for the care team**: Neutral prompts for clarification or shared decision-making.
6. **Safety and escalation note**: Any urgent review need, with no unsupported emergency claim.
7. **Quality review record**: Material facts checked, simplification risks found, comparison differences, and unresolved discrepancies.
8. **Professional review status**: `Required`, `Completed`, or `Not applicable`, with a brief rationale.

Do not present the output as a diagnosis, treatment plan, substitute for clinical care, or assurance that a finding is harmless. Keep provenance clear by labeling source-derived content, explanatory context, and reviewer-added changes.

## References

- PubMed: [Comparing large language models and human experts in interpreting MRI reports for personalized patient education](https://pubmed.ncbi.nlm.nih.gov/41865475/)
- PubMed: [Artificial Intelligence to Improve Patient-Physician Communication in Surgery](https://pubmed.ncbi.nlm.nih.gov/41948164/)
