# Tutorial: Anthropic Healthcare Coworker Workflows

**Stack:** Anthropic Claude for Healthcare (January 2026)
**Focus:** Complex Reasoning & Auditability
**Key Feature:** XML Tags & Thinking Blocks

## Introduction

Anthropic's Claude for Healthcare models are designed to be "Coworkers" - they excel at:
1. **Thinking Out Loud:** generating `<thinking>` blocks before final answers
2. **XML Protocol:** Using strict XML tags for structured outputs
3. **Evidence Synthesis:** Combining clinical guidelines, policies, and patient data
4. **Audit Trails:** Every decision is traceable and reviewable

## Available Coworkers

### 1. FHIR Development Coworker

Generate FHIR R4 compliant resources for healthcare interoperability:

```python
from fhir_development_coworker import FHIRDevelopmentCoworker

coworker = FHIRDevelopmentCoworker()

# Create a Patient resource
patient = coworker.generate_resource("Patient", {
    "identifier": "MRN-12345",
    "family_name": "Smith",
    "given_names": ["John"],
    "gender": "male",
    "birth_date": "1985-03-15",
})

print(patient["resource"])  # FHIR R4 Patient JSON
print(patient["trace"])     # Claude-style reasoning trace
```

### 2. Claims Appeals Coworker

Process insurance claim denials with evidence synthesis:

```python
from claims_appeals_coworker import ClaimsAppealsCoworker

coworker = ClaimsAppealsCoworker()

result = coworker.process_appeal(
    claim={"claim_id": "CLM-001", "procedure_code": "MRI-LUMBAR"},
    denial_reason="medical_necessity",
    clinical_records=[
        {"type": "clinical_note", "findings": ["Failed conservative therapy"]}
    ]
)

print(result["recommendation"])  # STRONG_APPEAL, MODERATE_APPEAL, or WEAK_APPEAL
print(result["appeal_letter"])   # Generated appeal letter
```

### 3. Care Coordination Coworker

Triage patient portal messages:

```python
from care_coordination_coworker import CareCoordinationCoworker

coworker = CareCoordinationCoworker()

result = coworker.triage_message({
    "message_id": "MSG-001",
    "patient_id": "PT-12345",
    "subject": "Chest pain",
    "content": "I've been having chest pain since yesterday",
})

print(result["urgency_level"])    # emergent, urgent, routine, low
print(result["routing"])          # Where to route the message
```

### 4. Research Literature Coworker

Search and synthesize biomedical literature from PubMed:

```python
from research_literature_coworker import ResearchLiteratureCoworker

coworker = ResearchLiteratureCoworker()

result = coworker.search_literature(
    query="SGLT2 inhibitors heart failure",
    filters={"year_min": 2024}
)

print(result["evidence_summary"])  # Evidence level distribution
print(result["synthesis"])         # Synthesized findings
```

### 5. Lab Results Coworker

Interpret lab results with patient-friendly explanations:

```python
from lab_results_coworker import LabResultsCoworker

coworker = LabResultsCoworker()

result = coworker.interpret_results([
    {"test_name": "Hemoglobin", "value": 10.5},
    {"test_name": "HbA1c", "value": 7.2},
])

print(result["patient_summary"])   # Patient-friendly summary
print(result["patterns"])          # Clinical patterns identified
```

### 6. Regulatory Drafter (Original)

The original regulatory drafting coworker:

```python
# In regulatory_drafter.py
# Simulates a "Co-worker" workflow with thinking phase
```

**The "Thinking" Phase:**
```xml
<thinking>
1. Analyze the Regulation (21 CFR).
2. Check the demographics in the clinical data.
3. Identify the discrepancy.
4. Decide to request a Waiver.
</thinking>
```

## Coworker Trace Format

All coworkers emit structured traces for auditing:

```xml
<thinking>Initial analysis of the input...</thinking>
<analysis>Detailed evaluation against criteria...</analysis>
<decision>Final determination with rationale...</decision>
```

## Integration with BioKernel

Register coworkers as MCP tools:

```python
from platform.biokernel import register_tool

register_tool("anthropic.fhir", FHIRDevelopmentCoworker())
register_tool("anthropic.claims_appeal", ClaimsAppealsCoworker())
register_tool("anthropic.care_triage", CareCoordinationCoworker())
register_tool("anthropic.literature", ResearchLiteratureCoworker())
register_tool("anthropic.lab_interpret", LabResultsCoworker())
```

## Running the Demos

```bash
# Run each coworker demo
python3 Skills/Anthropic_Health_Stack/FHIR_Development/coworker.py
python3 Skills/Anthropic_Health_Stack/Claims_Appeals/coworker.py
python3 Skills/Anthropic_Health_Stack/Care_Coordination/coworker.py
python3 Skills/Anthropic_Health_Stack/Research_Literature/coworker.py
python3 Skills/Anthropic_Health_Stack/Lab_Results/coworker.py
```

## Safety Considerations

Per Anthropic's acceptable use policy:
> "A qualified professional must review the content or decision prior to dissemination or finalization when Claude is used for healthcare decisions."

All coworkers flag critical values and support human-in-the-loop escalation.

---
*See [Anthropic_Health_STACK.md](../Anthropic_Health_STACK.md) for full documentation.*
