<!--
# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

-->

# Clinical Hallucination Detection Agent

**Version:** 1.0.0
**Status:** Production
**Date:** January 2026

## Overview

This agent implements continuous hallucination detection for clinical AI applications using the CHECK framework. It integrates structured clinical databases with information-theoretic classifiers to detect both factual and reasoning-based hallucinations in medical LLM outputs.

## Key Features

- **Factual Hallucination Detection**: Verifies clinical facts against trusted medical databases
- **Reasoning Hallucination Detection**: Identifies logical inconsistencies in medical reasoning
- **Continuous Learning**: Adapts to new clinical knowledge as databases update
- **Confidence Scoring**: Provides calibrated confidence scores for all predictions
- **Integration Ready**: Works with MedPrompt, clinical RAG, and other clinical AI pipelines

## Performance Benchmarks

Based on CHECK framework evaluation (arXiv:2506.11129):

| Metric | Value |
|--------|-------|
| Hallucination Detection Rate | 31% → 0.3% (with LLaMA3.3-70B) |
| AUC on MedQA (USMLE) | 0.95 |
| AUC on HealthBench | 0.96 |
| GPT-4o USMLE Improvement | +5 points (92.1%) |

## Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                  Hallucination Detection Pipeline            │
├─────────────────────────────────────────────────────────────┤
│                                                              │
│  ┌──────────────┐    ┌──────────────┐    ┌──────────────┐  │
│  │   LLM        │───→│   Claim      │───→│   Fact       │  │
│  │   Output     │    │   Extractor  │    │   Checker    │  │
│  └──────────────┘    └──────────────┘    └──────────────┘  │
│                                                   │         │
│                                                   ▼         │
│  ┌──────────────┐    ┌──────────────┐    ┌──────────────┐  │
│  │   Medical    │◀───│   Database   │◀───│   Reasoning  │  │
│  │   KB         │    │   Grounding  │    │   Validator  │  │
│  └──────────────┘    └──────────────┘    └──────────────┘  │
│                                                   │         │
│                                                   ▼         │
│                              ┌──────────────────────┐       │
│                              │   Confidence Score   │       │
│                              │   & Recommendations  │       │
│                              └──────────────────────┘       │
│                                                              │
└─────────────────────────────────────────────────────────────┘
```

## Components

### 1. Claim Extractor
Identifies verifiable medical claims from LLM outputs:
- Diagnosis statements
- Treatment recommendations
- Drug information
- Laboratory value interpretations
- Clinical guidelines citations

### 2. Fact Checker
Verifies claims against trusted sources:
- Medical literature (PubMed, Cochrane)
- Drug databases (DrugBank, FDA labels)
- Clinical guidelines (UpToDate, NICE)
- Medical ontologies (SNOMED-CT, ICD-10)

### 3. Reasoning Validator
Detects logical inconsistencies:
- Temporal reasoning errors
- Contraindication violations
- Dose calculation errors
- Differential diagnosis logic

### 4. Confidence Scorer
Provides calibrated probability estimates:
- Per-claim confidence scores
- Aggregate document confidence
- Uncertainty quantification

## Usage

```python
from hallucination_detector import HallucinationDetector, ClinicalKnowledgeBase

# Initialize with clinical knowledge base
kb = ClinicalKnowledgeBase()
detector = HallucinationDetector(knowledge_base=kb)

# Analyze LLM output
llm_response = """
Based on the patient's presentation, I recommend starting metformin 1000mg
twice daily for Type 2 diabetes management. This is the first-line treatment
per ADA guidelines 2024.
"""

# Run detection
result = detector.analyze(
    response=llm_response,
    source_context="65yo M with HbA1c 7.8%, no renal impairment",
    detection_mode="comprehensive"
)

# Review results
print(f"Confidence Score: {result.confidence:.2f}")
print(f"Hallucination Risk: {result.risk_level}")
print(f"Flagged Claims: {len(result.flagged_claims)}")

for claim in result.flagged_claims:
    print(f"  - {claim.text}")
    print(f"    Issue: {claim.issue_type}")
    print(f"    Confidence: {claim.confidence:.2f}")
```

## Configuration

```yaml
hallucination_detection:
  # Detection sensitivity
  sensitivity: "high"  # low, medium, high

  # Knowledge base sources
  knowledge_bases:
    - "pubmed"
    - "drugbank"
    - "uptodate"
    - "snomed_ct"

  # Confidence thresholds
  thresholds:
    low_risk: 0.9
    medium_risk: 0.7
    high_risk: 0.5

  # Reasoning validation
  reasoning:
    check_temporal: true
    check_contraindications: true
    check_dosing: true
    check_guidelines: true
```

## Integration with MedPrompt

```python
from medprompt_utils import MedPromptEngine
from hallucination_detector import HallucinationDetector

# Create pipeline
medprompt = MedPromptEngine()
detector = HallucinationDetector()

# Generate and validate
clinical_note = "..."
summary = medprompt.generate_clinical_summary(clinical_note)

# Check for hallucinations
validation = detector.analyze(summary, source_context=clinical_note)

if validation.risk_level == "high":
    # Regenerate with constraints
    summary = medprompt.generate_clinical_summary(
        clinical_note,
        constraints=validation.get_constraints()
    )
```

## Safety Considerations

1. **Human Review Required**: All high-risk flagged outputs must be reviewed by qualified clinicians
2. **Not Diagnostic**: This tool assists but does not replace clinical judgment
3. **Knowledge Currency**: Medical knowledge bases must be regularly updated
4. **Context Dependent**: Detection accuracy varies by clinical domain

## References

- [CHECK Framework](https://arxiv.org/html/2506.11129) - Continuous Hallucination Detection and Elimination
- [Nature Digital Medicine - Hallucination Framework](https://www.nature.com/articles/s41746-025-01670-7)
- [MedQA Benchmark](https://github.com/jind11/MedQA)
- [HealthBench Evaluation](https://healthbench.org)

## Dependencies

```
numpy>=1.24.0
sentence-transformers>=2.2.0
faiss-cpu>=1.7.0  # or faiss-gpu for acceleration
spacy>=3.5.0
scispacy>=0.5.0
transformers>=4.30.0
```

## License

MIT License - See repository root for full license.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->