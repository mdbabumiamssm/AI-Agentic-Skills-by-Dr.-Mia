# CAR-T Cell Design Agent

**ID:** `biomedical.immunology.car_t_design`
**Version:** 1.0.0
**Status:** Production
**Category:** Immunology / Cell Therapy / CAR-T

---

## Overview

The **CAR-T Cell Design Agent** assists in the rational design of chimeric antigen receptor (CAR) T cell therapies. CAR-T cells are engineered T lymphocytes expressing synthetic receptors that recognize tumor-associated antigens, representing a revolutionary approach to cancer immunotherapy with remarkable success in hematologic malignancies.

This agent integrates target antigen selection, scFv optimization, CAR architecture design, and safety feature engineering to support the development of next-generation CAR-T therapies with improved efficacy and reduced toxicity.

---

## Key Capabilities

### 1. CAR Architecture Components

| Domain | Function | Options |
|--------|----------|---------|
| **scFv** | Antigen recognition | Humanized, fully human, nanobody |
| **Hinge** | Flexibility | CD8α, CD28, IgG4 |
| **Transmembrane** | Membrane anchor | CD8α, CD28, 4-1BB |
| **Costimulatory** | Signal 2 | CD28, 4-1BB, OX40, ICOS |
| **Signaling** | Signal 1 | CD3ζ |

### 2. CAR Generations

| Generation | Architecture | Features | Examples |
|------------|--------------|----------|----------|
| **1st** | scFv-CD3ζ | Basic signaling | Early trials |
| **2nd** | scFv-CD28/4-1BB-CD3ζ | + Costimulation | Kymriah, Yescarta |
| **3rd** | scFv-CD28-4-1BB-CD3ζ | Dual costimulation | Research |
| **4th (Armored)** | CAR + cytokine/checkpoint | Enhanced function | TRUCKs |
| **5th** | CAR + IL-2Rβ signaling | STAT3/5 activation | Research |

### 3. Safety Features

| Feature | Mechanism | Application |
|---------|-----------|-------------|
| **Suicide genes** | iCasp9, HSV-TK | Emergency elimination |
| **ON-switch** | Small molecule control | Conditional activation |
| **Logic gates** | AND/OR/NOT gates | Tumor specificity |
| **Inhibitory CAR** | iCAR | Normal tissue protection |
| **Synthetic Notch** | Sequential activation | Enhanced selectivity |

### 4. FDA-Approved CAR-T Products

| Product | Target | Indication | Generation |
|---------|--------|------------|------------|
| **Kymriah** | CD19 | ALL, DLBCL | 2nd (4-1BB) |
| **Yescarta** | CD19 | DLBCL | 2nd (CD28) |
| **Tecartus** | CD19 | MCL | 2nd (CD28) |
| **Breyanzi** | CD19 | LBCL | 2nd (4-1BB) |
| **Abecma** | BCMA | Myeloma | 2nd (4-1BB) |
| **Carvykti** | BCMA | Myeloma | 2nd (4-1BB) |

---

## Usage

### Example Prompt

```text
Design a CAR-T therapy targeting HER2 for solid tumors with enhanced safety features.

Requirements:
- Target: HER2 (ERBB2) for breast cancer
- Safety concern: HER2 expression on normal tissues (heart, lung)
- Need: Logic-gated approach to improve tumor specificity
- Efficacy: Ability to persist and function in solid tumor TME

Tasks:
1. Select optimal anti-HER2 scFv
2. Design CAR architecture
3. Incorporate safety mechanisms
4. Add TME-resistance features
5. Generate complete construct sequence
6. Propose manufacturing and clinical testing approach
```

### Expected Output

```
## CAR-T Design Report: HER2-Targeted Solid Tumor Therapy

### Target Assessment

#### HER2 Expression Profile

| Tissue | Expression Level | Clinical Concern |
|--------|------------------|------------------|
| Breast cancer (HER2+) | Very high (3+) | Target |
| Gastric cancer | High (2-3+) | Target |
| Cardiac tissue | Low (1+) | Safety concern |
| Lung epithelium | Low-moderate | Safety concern |
| GI epithelium | Low | Minimal concern |

**Challenge:** On-target, off-tumor toxicity
**Strategy:** Logic-gated AND gate CAR system

### scFv Selection

#### Candidate Anti-HER2 scFvs

| scFv | Source | Kd (nM) | Epitope | Safety Profile |
|------|--------|---------|---------|----------------|
| 4D5 (trastuzumab) | Humanized | 0.1 | Domain IV | Cardiotoxicity risk |
| FRP5 | Humanized | 5.0 | Domain IV | Moderate affinity |
| **4D5-low** | Engineered | 50 | Domain IV | Reduced off-tumor |
| scFv-8 | Fully human | 8.2 | Domain II | Novel epitope |

**Selected:** 4D5-low affinity variant (Kd ~50 nM)
**Rationale:** Reduced binding to low-HER2 normal tissues while retaining tumor recognition

### CAR Architecture Design

#### Primary HER2-CAR (AND-Gate System)

```
┌─────────────────────────────────────────────────────────────┐
│                                                             │
│  RECEPTOR 1: HER2 Recognition (synNotch)                    │
│  ┌──────┐    ┌──────┐    ┌──────────────┐                   │
│  │ scFv │────│ Notch│────│ TF Release   │                   │
│  │4D5low│    │ core │    │ (Gal4-VP64)  │                   │
│  └──────┘    └──────┘    └──────┬───────┘                   │
│                                 │                           │
│                                 ▼                           │
│                    ┌────────────────────────┐               │
│                    │ UAS Promoter Activation│               │
│                    └────────────┬───────────┘               │
│                                 │                           │
│  RECEPTOR 2: Tumor-Specific CAR │                           │
│  ┌──────┐    ┌────┐  ┌────┐  ┌─▼───┐  ┌─────┐               │
│  │ scFv │────│ TM │──│4-1BB│──│CD3ζ │──│iCasp9│             │
│  │EpCAM │    │    │  │    │  │     │  │     │               │
│  └──────┘    └────┘  └────┘  └─────┘  └─────┘               │
│                                                             │
└─────────────────────────────────────────────────────────────┘

AND-Gate Logic:
- HER2 binding → synNotch activation → TF release
- TF activates CAR expression
- CAR recognizes EpCAM (tumor antigen 2)
- Both HER2 AND EpCAM required for T cell activation
```

#### Safety Features Incorporated

| Feature | Implementation | Function |
|---------|---------------|----------|
| **Low-affinity scFv** | 4D5-low (Kd 50nM) | Reduces off-tumor binding |
| **AND-gate logic** | synNotch + CAR | Dual antigen requirement |
| **iCasp9 suicide** | AP1903-inducible | Emergency elimination |
| **IL-15 knockdown** | Constitutive shRNA | Prevents uncontrolled expansion |
| **PD-1 knockout** | CRISPR edit | TME resistance |

### Complete Construct Design

#### Construct 1: synNotch Receptor

```
Component           Sequence Element                    Length
──────────────────────────────────────────────────────────────
Signal peptide      CD8α leader                         21 aa
scFv                4D5-low (VH-linker-VL)             245 aa
Hinge               CD8α                                45 aa
Notch core          Human Notch1 (1427-1752)           325 aa
Transcription factor Gal4-VP64                         216 aa
──────────────────────────────────────────────────────────────
Total                                                  852 aa
```

#### Construct 2: Inducible CAR

```
Component           Sequence Element                    Length
──────────────────────────────────────────────────────────────
Promoter            5x UAS - minimal CMV                 -
Signal peptide      GM-CSF leader                       22 aa
scFv                Anti-EpCAM (MOC31)                 247 aa
Hinge               IgG4-Fc mutant (L235E, N297Q)       45 aa
Transmembrane       CD28 TM                             27 aa
Costimulatory       4-1BB (137-169)                     33 aa
Signaling           CD3ζ (52-164)                      113 aa
──────────────────────────────────────────────────────────────
Total CAR                                              487 aa

Suicide gene        iCasp9                             296 aa
Linker              T2A                                 22 aa
──────────────────────────────────────────────────────────────
Total construct                                        805 aa
```

### TME Resistance Engineering

| Modification | Target | Mechanism |
|--------------|--------|-----------|
| **PD-1 KO** | PDCD1 | CRISPR editing |
| **TGF-β DN receptor** | TGF-β signaling | Dominant negative |
| **IL-15 tethered** | Autocrine support | Membrane-bound IL-15 |
| **Hypoxia-resistant** | HIF-1α stabilization | Constitutive activity |

### Manufacturing Considerations

#### Vector Design

| Element | Choice | Rationale |
|---------|--------|-----------|
| Vector type | Lentiviral (SIN) | Stable integration, safety |
| Promoter | EF1α | Strong, constitutive |
| Insert size | ~3.8 kb (dual construct) | Within packaging limit |
| Titer target | >10^8 TU/mL | Efficient transduction |

#### Process Overview

```
Leukapheresis
    ↓
T cell isolation (CD3+ selection)
    ↓
Activation (CD3/CD28 beads + IL-2)
    ↓
Transduction (MOI 5-10, Day 1)
    ↓
Expansion (IL-7/IL-15, 9-14 days)
    ↓
Harvest & QC
    ↓
Cryopreservation
    ↓
Release testing
    ↓
Infusion
```

#### Release Criteria

| Test | Specification |
|------|--------------|
| Viability | >70% |
| CAR+ cells | >20% |
| CD3+ purity | >95% |
| Sterility | Negative |
| Mycoplasma | Negative |
| Endotoxin | <5 EU/mL |
| Vector copy number | <5 copies/cell |
| Potency (cytotoxicity) | >30% lysis at 10:1 E:T |

### Preclinical Testing Plan

| Study | Model | Endpoints |
|-------|-------|-----------|
| In vitro cytotoxicity | HER2+/EpCAM+ cell lines | EC50, specificity |
| Logic gate validation | Mixed antigen cells | AND-gate function |
| In vivo efficacy | NSG mice, HER2+ xenograft | Tumor regression |
| Safety (on-target off-tumor) | HER2 transgenic mice | Organ toxicity |
| Suicide gene function | In vivo, AP1903 treatment | CAR-T depletion |

### Clinical Development Plan

#### Phase I Design (Dose Escalation)

| Dose Level | CAR+ T cells | Patients |
|------------|--------------|----------|
| DL1 | 1 x 10^7 | 3 |
| DL2 | 5 x 10^7 | 3 |
| DL3 | 1 x 10^8 | 3 |
| DL4 | 3 x 10^8 | 3-6 |

**Primary endpoints:** Safety, DLT, MTD
**Secondary endpoints:** ORR, PFS, CAR-T persistence
**Correlatives:** synNotch activation, CAR expression kinetics
```

### LLM Agent Integration

```python
@tool
def select_target_antigen(
    cancer_type: str,
    expression_database: str = "protein_atlas",
    safety_threshold: float = 2.0
) -> str:
    """
    Selects optimal CAR target antigen with safety assessment.

    Args:
        cancer_type: Target cancer type
        expression_database: Reference database
        safety_threshold: Tumor/normal expression ratio

    Returns:
        Ranked target antigens with safety profiles
    """
    pass


@tool
def design_scfv(
    target_antigen: str,
    affinity_range: tuple = (1, 100),
    format: str = "VH-linker-VL"
) -> str:
    """
    Designs/selects scFv for CAR targeting.

    Args:
        target_antigen: Target protein
        affinity_range: Desired Kd range (nM)
        format: scFv format

    Returns:
        scFv candidates with sequences and properties
    """
    pass


@tool
def design_car_architecture(
    scfv: str,
    generation: int = 2,
    costimulatory: str = "4-1BB",
    safety_features: list[str] = None
) -> str:
    """
    Designs complete CAR architecture.

    Args:
        scfv: scFv sequence or ID
        generation: CAR generation (2, 3, 4, 5)
        costimulatory: Costimulatory domain
        safety_features: Safety features to include

    Returns:
        Complete CAR design with sequences
    """
    pass


@tool
def design_logic_gated_car(
    antigen1: str,
    antigen2: str,
    logic: str = "AND"
) -> str:
    """
    Designs logic-gated CAR system.

    Args:
        antigen1: First target antigen
        antigen2: Second target antigen
        logic: Logic gate type (AND, OR, NOT)

    Returns:
        Logic-gated CAR system design
    """
    pass


@tool
def generate_manufacturing_protocol(
    car_construct: str,
    vector_type: str = "lentiviral",
    scale: str = "clinical"
) -> str:
    """
    Generates CAR-T manufacturing protocol.

    Args:
        car_construct: CAR construct design
        vector_type: Vector type (lentiviral, retroviral, non-viral)
        scale: Manufacturing scale

    Returns:
        Complete manufacturing protocol
    """
    pass
```

---

## Prerequisites

### Required Databases

| Database | Purpose |
|----------|---------|
| **Human Protein Atlas** | Expression profiling |
| **UniProt** | Protein sequences |
| **PDB** | Structural data |
| **IMGT** | Antibody sequences |

### Dependencies

```
biopython>=1.80
pandas>=2.0
numpy>=1.24
requests>=2.28
```

---

## Methodology

### CAR Design Workflow

```
Target Selection
├── Expression profiling (tumor vs normal)
├── Surface accessibility assessment
└── Safety evaluation
    ↓
scFv Selection/Engineering
├── Antibody screening
├── Affinity optimization
└── Humanization
    ↓
CAR Architecture Design
├── Hinge/TM selection
├── Costimulatory domain
└── Signaling domain
    ↓
Safety Engineering
├── Suicide genes
├── Logic gates
└── Regulatory elements
    ↓
Vector Construction
├── Promoter selection
├── Insert optimization
└── Production vector
    ↓
Manufacturing Protocol
├── T cell isolation
├── Transduction/editing
└── Expansion
    ↓
Release Testing
├── Phenotype
├── Potency
└── Safety
```

---

## Clinical Applications

### Hematologic Malignancies
- CD19+ B cell malignancies (ALL, DLBCL)
- BCMA+ multiple myeloma
- CD22, CD33 for AML

### Solid Tumors
- HER2+ breast/gastric cancer
- Mesothelin+ mesothelioma/ovarian
- GD2+ neuroblastoma
- EGFR+ glioblastoma

### Emerging Targets
- Claudin18.2 (gastric)
- CLDN6 (ovarian)
- GPC3 (HCC)
- ROR1 (various)

---

## Related Skills

- **TCR Engineering Agent:** TCR-based cell therapy
- **Neoantigen Vaccine Agent:** Target identification
- **Immune Repertoire Agent:** T cell characterization
- **Gene Editing Agent:** CRISPR modifications

---

## References

- **June et al. (2018):** "CAR T cell immunotherapy for human cancer." *Science*
- **Mackall et al. (2022):** "Next-generation CAR T cells: Overcoming challenges in solid tumors." *Nature Reviews Clinical Oncology*
- [CAR T Cell Clinical Trials](https://clinicaltrials.gov/)

---

## Author

**MD BABU MIA**
*Artificial Intelligence Group*
*Icahn School of Medicine at Mount Sinai*
md.babu.mia@mssm.edu
