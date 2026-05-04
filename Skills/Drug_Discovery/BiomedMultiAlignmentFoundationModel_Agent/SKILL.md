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
name: 'biomed-multi-alignment-foundation-model'
description: 'Use IBM biomed.omics.bl.sm.ma-ted-458m, a 458M-parameter foundation model trained on 2B+ biological samples across proteins, small molecules, and single-cell gene data, for cross-modal drug discovery tasks.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Biomed Multi-Alignment Foundation Model Agent

## Overview

This skill operationalizes IBM's `biomed.omics.bl.sm.ma-ted-458m` foundation model — a 458M-parameter multi-modality alignment (MaMMaL / MAMMAL) transformer trained on over 2 billion biological samples spanning proteins, small molecules, and single-cell gene expression. It exposes joint embeddings and task heads so an agent can answer cross-modal drug-discovery questions (e.g., target identification, ligand-binding prediction, transcriptomic perturbation) without stitching together separate single-modality models.

## When to Use This Skill

- Tasks that require a **shared embedding space** across proteins, ligands (SMILES), and gene-expression vectors.
- **Cross-modal retrieval**: find candidate small molecules for a target protein, or rank proteins likely to drive a transcriptomic phenotype.
- **Drug–target interaction (DTI)** scoring where both sequence and chemistry matter.
- **Virtual screening** seeded by a single-cell gene signature (e.g., disease vs. healthy state).
- **Modality transfer**: use a pretrained MAMMAL backbone to fine-tune a downstream head with limited labeled data.
- When existing single-modality skills (e.g., `Antibody_Design`, `Protein_Structure`, `scFoundation`) are insufficient because the question crosses modalities.

## Core Capabilities

1. **Model loading** — Pull `ibm/biomed.omics.bl.sm.ma-ted-458m` weights from Hugging Face and instantiate the MAMMAL encoder/decoder via the `biomed-multi-alignment` repo's tokenizer + task-prompt pipeline.
2. **Multi-modal tokenization** — Encode protein sequences (amino acids), small molecules (SMILES), and gene-expression profiles using the unified tokenizer the model was trained with, so inputs share one vocabulary.
3. **Cross-modal embedding extraction** — Produce dense vectors for any of the three modalities and project them into the shared latent space for similarity / retrieval.
4. **Task-prompt inference** — Use prompt-based task specification (the model is trained as a sequence-to-sequence learner over typed tokens) to run downstream tasks such as DTI, property prediction, or sequence labeling without retraining.
5. **Fine-tuning hooks** — Attach lightweight task heads (LoRA / adapter / linear probe) on top of frozen encoder representations for new datasets.
6. **Batch inference & caching** — Stream long candidate libraries (e.g., ChEMBL, UniProt subsets) through the model with embedding caches to keep cost bounded.
7. **Repository-backed Ma-TED workflow** — Use the `BiomedSciAI/biomed-multi-alignment` repository as the setup reference for `ibm/biomed.omics.bl.sm.ma-ted-458m`; expect notebook-oriented examples and dependency setup before extracting cross-modal embeddings for proteins, small molecules, and single-cell gene data in drug-discovery workflows such as target-ligand retrieval, molecule ranking against biological context, and single-cell signature alignment.

## Inputs / Outputs

**Inputs**
- Protein sequences (FASTA or raw amino-acid strings).
- Small molecules as SMILES.
- Single-cell or bulk gene-expression vectors (gene → expression value), aligned to the model's gene vocabulary.
- Task prompt (e.g., DTI scoring, property regression, retrieval).
- Optional: labeled dataset for fine-tuning, target list for screening.

**Outputs**
- Per-input embeddings in the shared latent space (numpy / torch tensors).
- Task-specific predictions (binding score, property value, ranked candidate list, cross-modal nearest neighbors).
- Optional fine-tuned checkpoint plus evaluation metrics on a held-out split.
- A run report summarizing inputs, model version, prompt template, and top-k results.

## References

- Source repository (finding): https://github.com/BiomedSciAI/biomed-multi-alignment
- Model card on Hugging Face: https://huggingface.co/ibm/biomed.omics.bl.sm.ma-ted-458m
- IBM Research project page (BiomedSciAI org): https://github.com/BiomedSciAI
- Related single-modality foundation work for context:
  - ESM-2 (protein language model): https://github.com/facebookresearch/esm
  - scFoundation (single-cell foundation model): https://github.com/biomap-research/scFoundation
  - MoLFormer (small-molecule foundation model): https://github.com/IBM/molformer
