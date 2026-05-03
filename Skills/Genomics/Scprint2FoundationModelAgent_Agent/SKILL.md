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
name: 'scprint2-foundation-model-agent'
description: 'Agent skill for applying the scPRINT-2 next-generation single-cell foundation model from the Cantini Lab to embedding, annotation, and inference tasks on scRNA-seq data.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# scPRINT-2 Single-Cell Foundation Model Agent

## Overview

This skill operationalizes the **scPRINT-2** single-cell foundation model from the Cantini Lab as an agentic capability for downstream single-cell analysis. It guides setup, model loading, inference, and integration of scPRINT-2 outputs (cell embeddings, annotations, gene-network signals) into broader genomics pipelines. Use it when an analyst needs a next-generation single-cell foundation model that complements (rather than replaces) the existing `scFoundation_Model_Agent` skill in this collection.

## When to Use This Skill

- The user asks for a "single-cell foundation model" and explicitly wants scPRINT-2 (or "scPRINT") rather than scFoundation, Geneformer, or scGPT.
- A workflow involves embedding scRNA-seq cells into a foundation-model latent space for clustering, label transfer, or batch integration.
- The user wants to perform zero-shot or fine-tuned cell type annotation using a pretrained transformer-style single-cell model.
- A pipeline requires gene-network or perturbation-aware inference signals on top of an AnnData object.
- Benchmarking scPRINT-2 against other single-cell foundation models on a held-out dataset.
- Reproducing or extending experiments published by the Cantini Lab around scPRINT-2.

## Core Capabilities

1. **Environment provisioning** — Provide reproducible setup steps for cloning `cantinilab/scPRINT-2`, creating a Python environment, and installing required dependencies (PyTorch, AnnData/Scanpy, Lightning, plus any model-specific extras documented in the upstream README).
2. **Checkpoint acquisition** — Retrieve the official scPRINT-2 model weights as published by the Cantini Lab (Hugging Face or repository release artifacts) and verify that the checkpoint loads successfully before any downstream call.
3. **Data preparation** — Validate that the input is a properly formed AnnData (`.h5ad`) object with raw counts, gene symbols matching the model's vocabulary, and required `obs` metadata; perform gene-symbol harmonization where needed.
4. **Cell embedding** — Run scPRINT-2 in inference mode to produce per-cell latent embeddings suitable for UMAP/PCA visualization, clustering (Leiden), and downstream classifiers.
5. **Cell-type annotation** — Apply scPRINT-2's prediction heads (or a lightweight fine-tuned classifier on its embeddings) to assign cell-type labels, returning calibrated confidence scores.
6. **Gene / network inference** — Where supported by the released model, extract gene-program or gene-network signals (e.g., attention-derived or decoder-derived) for use in regulatory analyses.
7. **Integration with downstream skills** — Hand off embeddings and annotations to other skills in this collection (visualization, differential expression, trajectory inference) via standard AnnData fields (`obsm`, `obs`).
8. **Benchmark harness** — Provide a thin wrapper for evaluating scPRINT-2 against `scFoundation_Model_Agent` and other single-cell foundation models on the same input dataset, reporting standard metrics (ARI, NMI, kBET, ASW, label-transfer accuracy).

## Inputs / Outputs

**Inputs**
- An `AnnData` object (`.h5ad`) with raw or lightly normalized scRNA-seq counts, valid gene identifiers, and optional ground-truth labels in `obs`.
- A model checkpoint identifier or local path for scPRINT-2 (release tag or Hugging Face repo).
- Optional configuration: target task (`embed`, `annotate`, `network`), batch size, device (`cuda` / `cpu`), and a results output directory.
- Optional reference label set or marker-gene panel for evaluation / annotation.

**Outputs**
- Updated `AnnData` with:
  - `adata.obsm["X_scprint2"]` — per-cell embedding matrix.
  - `adata.obs["scprint2_celltype"]` and `adata.obs["scprint2_celltype_score"]` — predicted labels and confidences (when annotation is run).
  - Optional `adata.uns["scprint2_network"]` — gene-program or gene-network outputs.
- A run manifest (JSON) recording model version, checkpoint hash, dataset shape, runtime, and key hyperparameters.
- Optional benchmark report comparing scPRINT-2 to alternative foundation models on the same input.

## References

- Source repository: <https://github.com/cantinilab/scPRINT-2>
- Cantini Lab (organization page): <https://github.com/cantinilab>
- Related skill in this collection: `Skills/Genomics/scFoundation_Model_Agent` for cross-model comparison.
- AnnData / Scanpy ecosystem (input format and downstream tooling): <https://github.com/scverse/scanpy>
