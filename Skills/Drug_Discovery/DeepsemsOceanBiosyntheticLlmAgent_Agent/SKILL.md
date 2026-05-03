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
name: 'deepsems-ocean-biosynthetic-llm-agent'
description: 'Agent that applies the DeepSeMS large language model to mine biosynthetic gene clusters and secondary metabolite potential from global ocean microbiome metagenomes.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# DeepSeMS Marine Microbiome Biosynthetic LLM Agent

## Overview
This skill operationalizes DeepSeMS, a large language model approach introduced in *Nature Computational Science* (2026) for uncovering the hidden biosynthetic potential of the global ocean microbiome. It guides an agent through ingesting marine metagenomic assemblies, running LLM-based prediction of biosynthetic gene clusters (BGCs) and secondary metabolite classes, and ranking novel natural product candidates for downstream chemistry follow-up. The goal is to accelerate marine drug discovery by surfacing previously unannotated BGCs that conventional rule-based tools (antiSMASH, DeepBGC) miss.

## When to Use This Skill
- You have metagenomic contigs, MAGs, or assembled scaffolds from ocean samples (e.g., Tara Oceans, Bio-GO-SHIP) and want to predict BGCs.
- Conventional BGC predictors return low yield or miss novel cluster classes in marine microbiomes.
- You need to prioritize candidate clusters for heterologous expression or molecular networking.
- You are building a marine natural products pipeline and want LLM-based annotation alongside antiSMASH / GECCO.
- You want to compare biosynthetic richness across ocean depth, latitude, or temperature gradients.
- You are exploring uncultured taxa where reference-based annotation is sparse.

## Core Capabilities
1. **Metagenome ingestion and preprocessing** — Accepts FASTA contigs/MAGs, performs ORF calling (Prodigal/Pyrodigal), generates protein sequences, and chunks them for LLM context windows while preserving genomic adjacency.
2. **DeepSeMS BGC inference** — Invokes the DeepSeMS model to predict BGC boundaries, core biosynthetic gene roles, and metabolite class probabilities (NRPS, PKS, RiPP, terpene, hybrid) directly from protein/token sequences.
3. **Novelty scoring against known clusters** — Compares predicted BGCs against MIBiG and antiSMASH reference clusters using sequence and embedding similarity to flag putatively novel chemistry.
4. **Cross-tool consensus** — Optionally cross-checks DeepSeMS predictions with antiSMASH 7, GECCO, and DeepBGC outputs and produces an ensemble confidence score.
5. **Ecological annotation** — Joins predicted BGCs with sample metadata (depth, temperature, oxygen, location) to map biosynthetic hot spots across ocean provinces.
6. **Candidate prioritization** — Ranks BGCs by novelty, completeness, host taxonomy, and predicted bioactivity class, emitting a shortlist suitable for synthetic biology follow-up.
7. **Reproducible reporting** — Produces per-cluster GenBank/JSON records, summary tables, and a provenance log capturing model version, thresholds, and inputs.

## Inputs / Outputs

### Inputs
- Assembled metagenomic FASTA files (contigs or MAGs) from marine samples.
- Optional sample metadata (sampling site, depth, environmental variables) as TSV/CSV.
- Optional reference databases: MIBiG, antiSMASH DB, NCBI taxonomy.
- Configuration: minimum contig length, BGC probability threshold, novelty cutoff, target metabolite classes.

### Outputs
- Per-sample table of predicted BGCs with coordinates, predicted class, probability, novelty score, and host MAG.
- GenBank or JSON files for each predicted BGC suitable for antiSMASH visualization or synthetic biology pipelines.
- Ranked candidate list for downstream cloning / heterologous expression.
- Ecological summary linking BGC density and class distribution to environmental gradients.
- Run manifest with model version, parameters, and input checksums.

## References
- Source paper: Xu T, Yang Y, Zhu R, Lin W, Li J. *DeepSeMS: revealing the hidden biosynthetic potential of the global ocean microbiome with a large language model.* Nat Comput Sci. 2026 Apr 30. https://pubmed.ncbi.nlm.nih.gov/42062603/
- antiSMASH 7 (rule-based BGC annotation): https://antismash.secondarymetabolites.org/
- MIBiG repository of characterized BGCs: https://mibig.secondarymetabolites.org/
- DeepBGC (deep-learning BGC detection): https://github.com/Merck/deepbgc
- GECCO (Gene Cluster prediction with Conditional Random Fields): https://github.com/zellerlab/GECCO
- Tara Oceans expedition data portal: https://fondationtaraocean.org/en/expedition/tara-oceans/
