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
name: bio-tcr-bcr-analysis-vdjtools-analysis
description: Calculate immune repertoire diversity metrics, compare samples, and track clonal dynamics using VDJtools. Use when analyzing repertoire diversity, finding shared clonotypes, or comparing immune profiles between conditions.
tool_type: cli
primary_tool: VDJtools
measurable_outcome: Execute skill workflow successfully with valid output within 15 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---

# VDJtools Analysis

## Basic Usage

```bash
# VDJtools requires Java
java -jar vdjtools.jar <command> [options]

# Or with wrapper script
vdjtools <command> [options]
```

## Calculate Diversity Metrics

```bash
# Basic diversity (Shannon, Simpson, Chao1, etc.)
vdjtools CalcDiversityStats \
    -m metadata.txt \
    output_dir/

# Metadata format (tab-separated):
# #file.name    sample.id    condition
# sample1.txt   S1           control
# sample2.txt   S2           treated
```

## Diversity Metrics Explained

| Metric | Description | Interpretation |
|--------|-------------|----------------|
| Shannon | Entropy-based diversity | Higher = more diverse |
| Simpson | Probability two random clones differ | 0-1, higher = diverse |
| InverseSimpson | 1/Simpson | Effective number of clones |
| Chao1 | Richness estimator | Total estimated clonotypes |
| Gini | Inequality coefficient | 0=equal, 1=dominated by one |
| d50 | Clones comprising 50% of repertoire | Lower = more oligoclonal |

## Sample Comparison

```bash
# Find overlapping clonotypes
vdjtools OverlapPair \
    -p sample1.txt sample2.txt \
    output_dir/

# Calculate overlap for all pairs
vdjtools CalcPairwiseDistances \
    -m metadata.txt \
    -i aa \
    output_dir/

# Overlap metrics: F2 (frequency-weighted Jaccard), Jaccard, MorisitaHorn
```

## Spectratype Analysis

```bash
# CDR3 length distribution (spectratype)
vdjtools CalcSpectratype \
    -m metadata.txt \
    output_dir/

# V/J gene usage
vdjtools CalcSegmentUsage \
    -m metadata.txt \
    output_dir/
```

## Clonal Tracking

```bash
# Track clones across timepoints
vdjtools TrackClonotypes \
    -m metadata_timecourse.txt \
    -x time \
    output_dir/

# Identify public clones (shared across individuals)
vdjtools JoinSamples \
    -m metadata.txt \
    -p \
    output_dir/
```

## Input Format

VDJtools accepts MiXCR output or standard format:

```
# Required columns (tab-separated):
count   frequency   CDR3nt  CDR3aa  V   D   J

# Example:
1500    0.15    TGTGCCAGC...    CASSF...    TRBV5-1*01  TRBD2*01    TRBJ2-7*01
```

## Convert from MiXCR

```bash
# Convert MiXCR output to VDJtools format
vdjtools Convert \
    -S mixcr \
    mixcr_clones.txt \
    output.txt
```

## Parse VDJtools Output in Python

```python
import pandas as pd

def load_diversity_stats(filepath):
    '''Load VDJtools diversity statistics'''
    df = pd.read_csv(filepath, sep='\t')
    return df

def load_overlap_matrix(filepath):
    '''Load pairwise overlap matrix'''
    df = pd.read_csv(filepath, sep='\t', index_col=0)
    return df

# Plot diversity across samples
def plot_diversity(stats_df, metric='shannon_wiener_index_mean'):
    import matplotlib.pyplot as plt

    plt.figure(figsize=(10, 6))
    plt.bar(stats_df['sample_id'], stats_df[metric])
    plt.xlabel('Sample')
    plt.ylabel(metric)
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig('diversity_plot.png')
```

## Related Skills

- mixcr-analysis - Generate input clonotype tables
- repertoire-visualization - Visualize VDJtools output
- immcantation-analysis - BCR-specific phylogenetics


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->