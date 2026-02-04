# 3D Genomics Skills

Skills for chromosome conformation capture (Hi-C) and 3D genome analysis.

## Sub-categories

| Category | Description | Skills |
|----------|-------------|--------|
| **hi-c-analysis** | Hi-C processing | Contact matrices, compartments, TADs, loops |

## Key Tools

- **HiC-Pro** - Hi-C data processing
- **cooler** - Contact matrix storage
- **cooltools** - Feature analysis
- **FAN-C** - Comprehensive Hi-C toolkit
- **HiCExplorer** - Visualization

## Capabilities

- Contact matrix generation and normalization
- A/B compartment calling
- TAD boundary detection
- Chromatin loop identification
- Differential interaction analysis
- 3D structure modeling

## Example Analysis

```python
import cooler
import cooltools

# Load Hi-C data
clr = cooler.Cooler('sample.mcool::resolutions/10000')

# Calculate expected
expected = cooltools.expected_cis(clr)

# Call compartments
eigenvector = cooltools.eigdecomp.eigs_cis(clr, phasing_track=gc_track)
```

---
*Source: bioSkills collection - integrated into BioMedical Skills Library*
