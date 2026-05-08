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
name: 'chromfound-scatac-foundation-model'
description: 'Apply ChromFound, a genome-wide foundation model for single-cell chromatin accessibility, to scATAC embedding, annotation, regulatory discovery, and validation.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# ChromFound scATAC Foundation Model

## Overview

This skill guides agents in applying ChromFound, a genome-wide foundation model for single-cell chromatin accessibility data, to scATAC-seq representation learning and downstream biological analysis. Use it to produce reproducible embeddings, transfer-learning plans, cell type annotations, candidate regulatory element hypotheses, and validation checks against conventional ArchR or Signac workflows.

The skill is designed for pragmatic analysis: inspect the available data and model assets first, run only documented or locally available ChromFound entry points, and treat model-derived outputs as hypotheses until compared with established quality control, clustering, peak-gene, motif, and marker-accessibility evidence.

## When to Use This Skill

- A user asks to use ChromFound or a scATAC foundation model on single-cell chromatin accessibility data.
- A project needs genome-wide scATAC embeddings for clustering, visualization, batch comparison, or integration.
- A user wants to transfer pretrained chromatin accessibility representations to a new tissue, disease, perturbation, or species-compatible dataset.
- A user needs cell type annotation support from scATAC embeddings and marker accessibility profiles.
- A project requires regulatory element discovery, enhancer prioritization, motif enrichment, peak-gene linkage, or differential accessibility hypotheses.
- A user asks to compare foundation-model outputs against ArchR, Signac, or other conventional scATAC workflows.
- A reviewer or collaborator asks for validation, reproducibility notes, or limitations of a ChromFound-based analysis.

## Core Capabilities

1. **Repository and environment triage**  
   Locate the ChromFound repository, notebooks, model weights, package requirements, and documented entry points before running analysis. Record missing dependencies, unavailable checkpoints, and assumptions rather than improvising unsupported APIs.

2. **Input inspection and quality control**  
   Verify accepted input formats such as fragment files, peak-by-cell matrices, AnnData, Seurat/Signac objects, Arrow files, or exported sparse matrices. Check basic scATAC quality metrics including cell counts, feature counts, sparsity, fragment depth, TSS enrichment when available, nucleosome signal when available, and batch labels.

3. **Embedding generation**  
   Run ChromFound inference or notebook workflows to generate cell-level and, when supported, region-level representations. Preserve model version, checkpoint path, genome build, feature space, preprocessing choices, and random seeds in the output notes.

4. **Cell state analysis and annotation**  
   Use embeddings for dimensionality reduction, neighborhood graph construction, clustering, and label transfer. Cross-check proposed labels with marker gene activity, marker peak accessibility, known lineage markers, and metadata before reporting annotations as final.

5. **Regulatory element discovery**  
   Prioritize candidate regulatory elements by combining model-derived representations with differential accessibility, motif enrichment, peak-to-gene links, enhancer annotations, and cell type specificity. Report candidates as ranked hypotheses with genomic coordinates, supporting evidence, and validation gaps.

6. **Transfer learning and adaptation**  
   Plan fine-tuning or feature extraction for new datasets by separating frozen embedding workflows from trainable adaptation workflows. Track train, validation, and held-out splits by donor, sample, batch, or biological condition to avoid leakage.

7. **Benchmarking against conventional workflows**  
   Compare ChromFound outputs with ArchR or Signac results for QC retention, low-dimensional structure, cluster stability, marker accessibility, differential accessibility, motif enrichments, and biological interpretability. Avoid claiming superiority unless the project supplies a controlled benchmark.

8. **Reproducible reporting**  
   Produce concise output artifacts: commands or notebooks used, environment details, input paths, output paths, summary metrics, plots generated, candidate regulatory elements, cell annotation tables, limitations, and next validation steps.

## Inputs / Outputs

**Inputs**

- ChromFound source repository, notebook, package, or local clone path.
- scATAC-seq data in a documented or convertible format, such as fragments, count matrices, AnnData, Seurat/Signac objects, ArchR Arrow projects, BED/peak files, or metadata tables.
- Genome build and feature definition, such as hg38, hg19, mm10, peaks, bins, genes, or model-specific region vocabulary.
- Optional pretrained checkpoint, tokenizer or region vocabulary, configuration file, and inference parameters.
- Optional reference labels, marker lists, ArchR outputs, Signac outputs, motif databases, gene annotations, or enhancer catalogs for validation.

**Outputs**

- ChromFound embedding files, reduced-dimension coordinates, cluster labels, and annotation tables when supported by the available workflow.
- Ranked candidate regulatory elements, cell type or condition associations, motif or peak-gene evidence, and differential accessibility summaries.
- Validation comparison against ArchR, Signac, or other baseline workflows, including agreement, discrepancies, and unresolved quality concerns.
- Reproducibility note containing repository URL or commit, commands or notebooks used, environment, input paths, output paths, model checkpoint, genome build, and known limitations.

## References

- ChromFound GitHub repository: https://github.com/SAIS-LifeScience/ChromFound
- ArchR GitHub repository: https://github.com/GreenleafLab/ArchR
- Signac GitHub repository: https://github.com/stuart-lab/signac
