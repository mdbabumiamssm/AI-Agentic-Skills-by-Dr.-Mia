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

---
name: single-cell-rna-qc
description: Perform quality control on single-cell RNA-seq data (.h5ad, 10x .h5, or 10x directories) using scverse best practices, MAD-based filtering (log1p counts/genes, high-tail MT%), and generate filtered AnnData plus QC plots and summary JSON. Use when users request scRNA-seq QC, filtering low-quality cells, data quality assessment, or QC visualizations.
measurable_outcome: Execute skill workflow successfully with valid output within 15 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---

# Single-Cell RNA-seq Quality Control

Automated QC workflow for single-cell RNA-seq data following scverse best practices.

## Approach 1: Complete QC Pipeline (Recommended)

Use the convenience script `scripts/qc_analysis.py` for end-to-end QC:

```bash
python3 scripts/qc_analysis.py input.h5ad
python3 scripts/qc_analysis.py raw_feature_bc_matrix.h5
python3 scripts/qc_analysis.py /path/to/10x_directory/
```

**When to use this approach:**
- Standard QC workflow with dataset-adaptive thresholds
- Batch processing or quick exploratory analysis
- Users want a reproducible, one-command pipeline

**Key parameters:**
- `--mad-counts`, `--mad-genes`, `--mad-mt` - MAD thresholds
- `--mt-threshold` - Hard mitochondrial % cutoff
- `--min-cells` - Gene filtering threshold
- `--mt-pattern`, `--ribo-pattern`, `--hb-pattern` - Gene patterns
- `--no-log1p` - Disable log1p transform for MAD on counts/genes

**Outputs (in `<input_basename>_qc_results/` by default):**
- `qc_metrics_before_filtering.png`
- `qc_filtering_thresholds.png`
- `qc_metrics_after_filtering.png`
- `<input_basename>_filtered.h5ad`
- `<input_basename>_with_qc.h5ad`
- `qc_summary.json`

## Approach 2: Modular Building Blocks (Custom Workflows)

For custom analysis workflows, use functions from `scripts/qc_core.py` and `scripts/qc_plotting.py`:

```python
import anndata as ad
from qc_core import calculate_qc_metrics, build_qc_masks, filter_cells

adata = ad.read_h5ad('input.h5ad')
calculate_qc_metrics(adata, inplace=True)

masks = build_qc_masks(
    adata,
    mad_counts=5,
    mad_genes=5,
    mad_mt=3,
    mt_threshold=8,
    counts_transform='log1p',
    genes_transform='log1p'
)

adata_filtered = filter_cells(adata, masks['pass_qc'])
```

**When to use this approach:**
- Non-standard filtering logic (subset-specific thresholds)
- Partial execution (metrics only, plots only)
- Integration with larger pipelines

## Best Practices

1. **Use log1p for counts/genes** - Stabilizes MAD thresholds for heavy-tailed distributions.
2. **Filter high MT% only** - Low MT% is usually not problematic.
3. **Inspect plots** - Validate that filtering aligns with biology and tissue context.
4. **Be permissive by default** - Preserve rare cell populations; filter further later if needed.

## Reference Materials

For deeper rationale, parameter guidance, and troubleshooting, see:
- `references/scverse_qc_guidelines.md`

## Next Steps After QC

- Ambient RNA correction (SoupX, CellBender)
- Doublet detection (scDblFinder, scrublet)
- Normalization and downstream analysis


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->