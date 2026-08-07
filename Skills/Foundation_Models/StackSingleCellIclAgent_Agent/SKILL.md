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
name: 'stack-single-cell-icl-agent'
description: 'Agent for Arc Institute Stack, a single-cell foundation model that performs in-context learning at inference time without per-task fine-tuning.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Stack In-Context Single-Cell Foundation Model Agent

## Overview

This skill operationalizes Arc Institute's **Stack**, a single-cell foundation model designed for **in-context learning (ICL)** at inference time. Unlike conventional single-cell foundation models that require task-specific fine-tuning, Stack accepts a small set of labeled "context" cells alongside query cells in a single forward pass and produces predictions directly. The agent helps users prepare AnnData inputs, assemble context/query batches, run inference, and interpret outputs for tasks such as cell-type annotation, perturbation response prediction, and condition transfer.

## When to Use This Skill

- The user wants to annotate or classify single-cell RNA-seq data using only a handful of labeled reference cells per class (few-shot setting).
- The user needs to predict cellular response to a new perturbation, condition, or tissue without retraining a model.
- The user is comparing in-context learning approaches against fine-tuned baselines such as scGPT, Geneformer, or scFoundation.
- The user is exploring transferable single-cell representations across datasets, batches, or species using prompt-style conditioning.
- The user is building an inference pipeline that must avoid the cost or operational complexity of per-task fine-tuning.

## Core Capabilities

1. **Environment bootstrap** — Guide cloning of `ArcInstitute/stack`, creation of an isolated Python environment, installation of dependencies, and verification of GPU/CUDA availability for transformer inference.
2. **Checkpoint acquisition** — Locate and download the released Stack model weights and tokenizer/gene vocabulary artifacts referenced in the repository, and validate file integrity before use.
3. **AnnData preparation** — Normalize raw counts, align gene symbols to the model's expected vocabulary, filter low-quality cells, and split data into context (labeled support) and query (unlabeled) sets.
4. **In-context prompt assembly** — Construct mixed batches of context and query cells with the labels, metadata, or perturbation tokens required by Stack's ICL interface, mirroring the patterns shown in the repository's Jupyter notebooks.
5. **Inference execution** — Run forward passes for cell-type prediction, perturbation response, or embedding extraction without weight updates, controlling batch size and context length to fit available GPU memory.
6. **Result post-processing** — Decode model outputs into per-cell predictions or embeddings, attach them back to the AnnData object, and export tables or UMAP-ready matrices for downstream analysis.
7. **Evaluation harness** — Compute few-shot accuracy, macro-F1, or perturbation correlation metrics by holding out labeled cells and varying the size of the in-context support set.
8. **Comparative benchmarking** — Optionally contrast Stack's zero-shot/few-shot performance with fine-tuned single-cell foundation model baselines on the same query split.
9. **Stack ICL workflow framing** — Treat each task as inference-only prompt conditioning rather than fine-tuning: package labeled reference/context cells and unlabeled query cells from AnnData together, keep benchmark comparisons matched on the same held-out query split and context budget, and validate predicted annotations against conventional annotation baselines before treating them as ground truth.

## Inputs / Outputs

**Inputs**
- An AnnData (`.h5ad`) file containing single-cell expression counts and, where available, cell-level labels or condition annotations.
- A small labeled support set (context cells) per class or condition of interest.
- A query set of unlabeled cells to be predicted or embedded.
- Stack model checkpoint(s) and gene vocabulary from the upstream repository.
- Optional configuration: context size, batch size, target task (annotation, perturbation, embedding), device.

**Outputs**
- Per-cell predictions (class labels, probabilities, or perturbation response vectors) written back into `adata.obs` / `adata.obsm`.
- Cell embedding matrices suitable for clustering, UMAP, or transfer to downstream classifiers.
- An evaluation report summarizing few-shot metrics across context-set sizes, when ground-truth labels are provided.
- Reproducible run logs capturing checkpoint hash, context composition, and inference parameters.

## References

- Stack repository (Arc Institute): https://github.com/ArcInstitute/stack
- Geneformer foundation model: https://github.com/jkobject/geneformer
- scGPT foundation model: https://github.com/bowang-lab/scGPT
- scFoundation foundation model: https://github.com/biomap-research/scFoundation
- AnnData / Scanpy ecosystem: https://github.com/scverse/scanpy
