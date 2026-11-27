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
name: 'single-cell-foundation-model-data-selection'
description: 'Compare pretraining dataset size and diversity for single-cell foundation models using scFM-dataselection notebooks, population separation, and zero-shot versus fine-tuned evaluation.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Single-Cell Foundation Model Data Selection

## Overview

Use the scFM-dataselection experiments to assess how pretraining subset size and sampling strategy affect downstream single-cell models. Build comparable subset manifests and notebook-backed comparisons of zero-shot and fine-tuned performance, retaining uncertainty and compute costs rather than assuming more cells improve every task. The repository's core comparison covers scVI, SSL, and Geneformer on scTab subsets, with classification and batch integration evaluations; verify additional model and task coverage against the selected revision.

## When to Use This Skill

- Choosing a pretraining corpus size under a compute or storage budget.
- Comparing uniform random sampling, cell-type reweighting, and geometric sketching.
- Testing whether broader population coverage improves transfer to a specific tissue or cohort.
- Comparing zero-shot classification with fine-tuned classification on consistent evaluation populations.
- Reproducing repository figures or auditing whether apparent scaling gains reflect leakage, changed splits, or unequal training budgets.

## Core Capabilities

1. **Establish a reproducible experiment contract.** Record the repository commit, dataset release, model checkpoint, preprocessing, task, sampling method, subset size, seeds, and available resources. Inspect actual script arguments and notebook input paths before execution: README examples include stale names and variable mismatches, including `random_downsampling.py` versus the listed `random_subsampling.py`. Resolve differences from code at the pinned revision, and record adaptations.

2. **Define comparable data subsets.** Follow `data/preprocess/` and `downsampling/` to select stable cell identifiers and materialize subsets without loading the full expression matrix unnecessarily. Compare sampling methods at equal realized training-cell counts; compare sizes within a method while holding architecture, gene representation, and evaluation population fixed. Use the repository's percentage and seed mapping when reproducing its experiments:

   | Corpus percentage | Replicate seeds |
   | --- | --- |
   | 1 | 0–4 |
   | 10 | 5–9 |
   | 25 | 10–14 |
   | 50 | 15–19 |
   | 75 | 20–24 |
   | 100 | 25–29 |

   Verify which combinations actually exist; do not treat every method at 100% as an independent corpus. Record realized cells, donors, studies, tissues, batches, and cell-type counts where metadata permits. Report rare-population coverage and abundance changes separately from dataset size. Geometric sketching uses an HVG/PCA preparation route; include its cost and distinguish expression-space coverage from biological diversity.

3. **Separate development and evaluation populations.** Audit overlap among pretraining, validation, fine-tuning, and final test sets using cell IDs and available donor/study identifiers. The repository's `downsampling/train_val_test_split.py` shuffles cell indices; it does not establish donor or study independence. For a new generalization experiment, define grouped holdouts before sampling and fit learned preprocessing on permitted training data. For exact reproduction, retain original splits and disclose their limitations instead of silently replacing them. Freeze downstream test populations across conditions, restrict tuning to development data, and label unknown pretraining overlap explicitly.

4. **Preserve model-specific input semantics.** Validate gene identifiers, order, missing genes, label mappings, and count transformations against the model loader and checkpoint. The preprocessing workflow produces raw H5AD, processed H5AD, and processed parquet representations; SSL evaluation requires matching the scTab gene space and zero-padding missing genes, while Geneformer requires its tokenizer and dictionaries. Compute fine-tuning class weights from training labels only. Do not feed a single normalized matrix indiscriminately to every architecture.

5. **Evaluate the repository's actual task paths.** Inspect `pretraining/`, `finetuning/`, and `eval/` for the chosen architecture and use their verified entry points. Compare zero-shot classification and fine-tuned classification separately, using the same downstream split and label policy; record any classifier fitting used by the zero-shot evaluator. Run zero-shot batch integration as a distinct task. Do not infer a fine-tuned integration experiment merely from the existence of a fine-tuning directory. Include supported PCA or untrained baselines with matching evaluation data, and record checkpoint provenance, optimizer steps, epochs, early stopping, and fine-tuning budget so data effects are not conflated with training exposure.

6. **Analyze results through the notebooks.** Start with `plotting/heatmaps.ipynb`, `plotting/lineplots.ipynb`, and `plotting/pretraining_dataset_sizes.ipynb`; inspect cells and required CSVs before running them. The heatmap notebook uses `micro_f1` for classification and `avg_bio` for zero-shot integration. Inspect evaluator definitions and retain batch-removal metrics alongside biological-conservation metrics rather than calling `avg_bio` an overall integration score. Use dataset-specific comparisons and replicate summaries; report missing conditions and avoid pooling incompatible metrics.

   The heatmap function selects the smallest nonzero percentage whose score exceeds 95% of the maximum observed score within a method/dataset. This is a descriptive rule, not evidence of statistical equivalence or a guaranteed plateau; inspect whether inputs are replicate-level or already aggregated. Report replicate dispersion and gains or regressions relative to the full-data reference. Check notebook variable bindings: the inspected fine-tuned spike-in classification block selects zero-shot dataframes despite its title. Correct that binding in an analysis copy before interpreting the figure and document the change.

7. **Quantify resource tradeoffs and select a defensible subset.** Report measured performance, uncertainty, population coverage, preprocessing time, training time, inference time, and storage for each viable condition. Prefer a smaller subset only when it meets a stated task-specific tolerance and coverage requirement; allow different choices for different tasks. Distinguish measured values from estimates and published outputs. Never manufacture missing scores or claim larger datasets, greater diversity, or a particular sampling method universally wins.

### Execution Scope and Compute Requirements

The 15-minute outcome applies to a bounded audit or analysis with prepared inputs, not full corpus download, installation, pretraining, and evaluation. Start by selecting one of these scopes:

- **Existing-results analysis:** With a working Python/Jupyter environment and local evaluation CSVs, validate the experiment manifest, inspect splits, and execute the relevant plotting/summary cells. CPU analysis is usually sufficient; measure local memory and runtime rather than promising a hardware minimum.
- **Checkpoint smoke evaluation:** With compatible data, dependencies, and a selected checkpoint already available, run a small evaluator batch to verify gene alignment, output shape, metric schema, and runtime. Mark its scores as smoke-test results, not the full experiment.
- **Full reproduction:** Prepare a separate resource estimate and run plan. The root README warns dependency installation alone can take more than 15 minutes. Its model collection is approximately 36.5 GB; download only required checkpoints when possible. The documented Geneformer pretraining example uses DeepSpeed on one node with eight GPUs; this is an example configuration, not a verified minimum for every model.

Check Python/package compatibility, model-specific forks, GPU/CUDA support, available GPU memory, CPU RAM, and disk capacity before training. Budget for raw and processed H5AD files, parquet copies, tokenized data, indices, checkpoints, and temporary outputs. Use backed/chunked reads where supported, and profile one representative batch before estimating full runtime; the inspected documentation does not establish universal RAM or GPU-memory requirements.

For inference comparisons, inspect `eval/timing_benchmark.py` and `plotting/inference_time.ipynb`. The latter reads `metrics_csvs/timing_benchmark_results.csv` and excludes the `Geneformer+tokenization` category. Report embedding-only and end-to-end timing separately if tokenization or preprocessing changes the comparison.

### Bounded Workflow and Acceptance Checks

1. During the initial inspection, pin source revisions, select a task/model pair, and inventory available manifests, checkpoints, CSVs, and notebook dependencies.
2. Validate subset sizes and uniqueness, gene order, population overlap, downstream split identity, and task/regime labels. Stop dependent evaluation if these checks fail; report the concrete mismatch.
3. Analyze existing results or execute the selected smoke evaluation within the remaining time. Use a clean notebook kernel and explicit paths; do not rely on stale displayed notebook output as newly reproduced evidence.
4. Produce a comparison table, notebook provenance, resource ledger, and conditional recommendation. If results or compute are unavailable, produce a validated experiment plan and explicit missing-input report; label the status `planned` or `blocked`, never `evaluated`.

A completed analysis has traceable input rows, correctly separated regimes, documented split checks, a baseline or explicit reason it is absent, and conclusions supported by available measurements. Record unverified donor/study overlap as a limitation. Distinguish exact repository reproduction from improved experimental controls or user-dataset extensions.

## Inputs / Outputs

**Inputs**

- Repository revision and selected notebook/evaluator paths.
- scTab or compatible expression data with gene metadata and stable cell identifiers; donor, study, tissue, batch, and cell-type annotations where available.
- Pretraining subset indices, downstream train/validation/test manifests, and any permitted population-overlap policy.
- Model checkpoints, preprocessing/tokenizer assets, evaluation CSVs, or an explicit request to generate them.
- Compute/storage limits, task priorities, performance tolerance, and requested execution scope.

**Outputs**

- Experiment manifest keyed by model, sampling method, percentage, realized training cells, seed, dataset, task, regime, and checkpoint revision.
- Split/coverage audit documenting duplicates, population overlap, excluded cells, missing metadata, and gene alignment.
- Notebook-backed comparison table containing metric definitions, values, replicate counts, variability, baseline differences, source CSV paths, and execution status.
- Resource ledger separating preprocessing, training, and inference costs, with measured versus estimated fields clearly marked.
- Conditional subset recommendation with supporting tradeoffs, unresolved limitations, and a reproducible next-run specification when evaluation is incomplete.

## References

- [scFM-dataselection source repository](https://github.com/microsoft/scFM-dataselection) — experiment overview, setup, checkpoints, and workflow entry points.
- [Data preprocessing instructions](https://github.com/microsoft/scFM-dataselection/blob/main/data/preprocess/README.md) — representations, percentage/seed mapping, and downstream inputs.
- [Downsampling workflow](https://github.com/microsoft/scFM-dataselection/tree/main/downsampling) — sampling implementations and split helper.
- [Pretraining instructions](https://github.com/microsoft/scFM-dataselection/blob/main/pretraining/README.md) — model preparation and example compute configuration.
- [Evaluation instructions](https://github.com/microsoft/scFM-dataselection/blob/main/eval/README.md) — classification, integration, and timing entry points; verify examples against code.
- [Analysis notebooks](https://github.com/microsoft/scFM-dataselection/tree/main/plotting) — figure logic and expected result schemas.
- [Evaluating the role of pre-training dataset size and diversity on single-cell foundation model performance](https://doi.org/10.1101/2024.12.13.628448) — associated preprint linked by the repository.
