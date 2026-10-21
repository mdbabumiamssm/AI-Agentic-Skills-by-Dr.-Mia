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
name: 'protein-geometry-agent-evaluation'
description: 'Evaluate protein-structure agents by separating metadata lookup from coordinate geometry, fixing structure selections and tool conditions, and reporting reproducible calculations with uncertainty-aware comparisons.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Protein Geometry Agent Evaluation

## Overview

Evaluate structural-biology agents using separate measures for recorded-fact lookup and calculations from protein coordinates. ProteinGeoBench reports that performance and apparent rankings depend on task family, repeated-run budget, harness, and access conditions; these findings motivate controlled evaluation rather than universal model rankings. This workflow produces auditable coordinate selections, independently checked answers, and comparisons that expose uncertainty.

## When to Use This Skill

- Compare agents answering questions about protein structures.
- Diagnose whether a wrong answer comes from retrieval, atom selection, geometry, or answer formatting.
- Evaluate changes to an agent harness, tool budget, or network policy.
- Reproduce a published structural-biology evaluation with verified task assets.
- Audit a claimed model advantage that may depend on task composition or too few runs.

## Core Capabilities

1. **Verify evidence and assets.** Separate published observations from local results, and establish the exact task, coordinate, and verifier versions before execution.
2. **Define structural selections.** Resolve structure identity, assembly, model, chains, residue numbering, alternate conformers, and atom membership explicitly.
3. **Reproduce geometry.** Compute quantities directly from selected coordinates with declared units, conventions, tolerances, and independent checks.
4. **Control evaluation conditions.** Freeze harness, tools, access, prompts, budgets, grading, and repeated-run policies so comparisons have interpretable scope.
5. **Quantify uncertainty.** Report task-family results and paired model differences without treating repeated trials as independent tasks.

### Workflow

#### 1. Scope a bounded evaluation

Use the 15-minute target for a verified pilot or analysis of existing rollouts, not a promise to finish an entire benchmark. Reserve approximately three minutes for asset checks, four for reference calculations, five for agent trials or log analysis, and three for validation and reporting.

Select a small, declared task subset containing both lookup and numeric work when possible. Freeze the selection before viewing model outcomes; record whether tasks were sampled, chosen for convenience, or selected for difficulty. If assets or execution time are insufficient, return a blocked or partial report with the exact missing prerequisites and completed coverage.

#### 2. Verify task assets before executable examples

Read the versioned paper and follow its benchmark link. The linked repository describes Harbor tasks and executable verifiers, but repository availability alone does not establish that an individual task is runnable.

Before proposing or running any benchmark command:

- Resolve and record the repository commit or release, task IDs, prompt contents, input files, dependency versions, and documented runner interface.
- Inspect each selected task's actual instructions, environment requirements, output contract, verifier, and reference-answer provenance.
- Confirm referenced structure files exist, parse successfully, contain the requested entities, and have recorded SHA-256 checksums.
- Inspect executable verifier code and dependencies before running it. Keep hidden answers and verifier internals outside the evaluated agent's accessible workspace.
- Check that the verifier accepts a known-correct answer and rejects a deliberately wrong value or malformed response.

Do not invent download paths, task IDs, CLI flags, expected answers, or claims of successful execution. If assets cannot be verified, provide a specification and verification checklist instead of purportedly runnable benchmark examples. Any synthetic fixture must be labeled as a local implementation check, not a ProteinGeoBench task.

#### 3. Freeze a task and selection manifest

Assign each task a primary family: metadata lookup, coordinate calculation, or mixed. Keep mixed tasks separate unless their components have predefined scoring rules. A numeric field copied from a database remains lookup; numeric formatting alone does not make a task geometric.

For every coordinate task, record:

| Field | Required decision |
| --- | --- |
| Structure | Accession or local identity, file format, revision or retrieval date, checksum |
| Coordinate scope | Asymmetric unit or named biological assembly; assembly transforms if applied |
| Model | Explicit model identifier; for ensembles, per-model calculation and aggregation rule |
| Chains | Exact identifiers and whether they are author or label identifiers; preserve mappings |
| Residues | Numbering namespace, residue name, sequence number, insertion code |
| Atoms | Atom names, element filters, hydrogens, ligands, waters, and heteroatom inclusion |
| Alternate locations | Task-prescribed conformer or deterministic, chemically coherent selection policy |
| Missing data | Missing residues or atoms, zero occupancy, and exclusion or failure policy |
| Geometry | Quantity, units, atom correspondence, cutoff, alignment, and periodicity conventions |

Never silently substitute the first chain, first model, nearby residue, or a homolog. Preserve multi-character chain IDs where present. Avoid combining incompatible alternate conformers; document occupancy-based selection and tie-breaking when used. Missing atoms are unavailable observations, not zero coordinates.

#### 4. Establish an independent geometric reference

Use a coordinate parser with a recorded version. Save the selected atom identifiers and coordinates or a reproducible selection expression plus their digest. Verify coordinate units before calculations; report distances in angstroms and angles in explicitly stated degrees or radians.

- **Distance:** For selected positions `a` and `b`, compute `sqrt(sum((a_i - b_i)^2))`. Specify whether a residue-level distance means C-alpha, centroid, or minimum interatomic distance.
- **Contacts:** Define atom filters, cutoff and strict/inclusive boundary, intra-chain exclusions, and whether counts represent unique atom pairs or residue pairs. Avoid self-pairs and double counting.
- **Angles and dihedrals:** Declare atom order and signed-angle convention. Reject coincident points or collinear configurations when the requested angle is undefined; compare periodic angles using wrapped differences.
- **RMSD:** Declare atom correspondence, matched count, and whether rigid-body superposition is performed. For aligned RMSD, use a proper rotation without reflection; report excluded atoms and alignment scope.
- **Clashes:** Declare radii source, overlap threshold, hydrogen treatment, and bonded-neighbor exclusions. A distance cutoff alone is not a universally defined clash criterion.

Check a small selected calculation through a second implementation or direct arithmetic. Check relevant invariants, such as distance symmetry and invariance under rigid translation and rotation. Do not replace coordinate-derived answers with remembered or retrieved literature values.

Record the unrounded result, output precision, and grading tolerance before inspecting agent answers. For ordinary scalar values, a declared rule may use `abs(predicted - reference) <= atol + rtol * abs(reference)`; justify tolerances using task semantics and coordinate precision. Preserve a verified benchmark's official grading rule, and report any supplementary geometric audit separately.

#### 5. Execute controlled agent trials

Freeze the model endpoint and date, agent harness and version, system/task prompts, sampling and reasoning settings, installed tools, dependency versions, token budget, wall-time limit, and tool-call limit. Record unsupported or unavailable settings rather than implying they were controlled.

Define each access condition explicitly:

- **Network enabled:** Record available services, credentials scope, and cached resources.
- **Domain restricted:** Record the exact blocked domains, subdomain handling, redirects, and enforcement mechanism; do not claim this removes all equivalent information sources.
- **Offline:** Disable network access at the execution boundary and document preloaded structures, packages, metadata, and caches.

Verify access conditions with benign connectivity checks inside the actual agent environment. Hold local assets constant across access contrasts unless asset availability is the factor being tested. Use fresh task workspaces and prevent answer leakage through prior transcripts or shared caches.

Run the same task set with the same planned number of independent attempts per condition. Record seeds where supported and every attempted run, including timeouts and tool failures. Randomize or interleave condition order when practical. Preserve stdout, stderr, tool traces, final answers, elapsed time, and usage.

Change one factor at a time for interpretable contrasts. A model's native CLI may also change prompts, defaults, and tools; label that comparison as a bundled configuration change unless those components are matched. Grade all compared conditions with the same frozen verifier or grader pass.

#### 6. Score and compare with uncertainty

Apply the predefined output parser and scoring rules. Classify errors as retrieval, wrong structure/selection, geometric calculation, units/convention, formatting, access failure, or execution failure; attach evidence rather than guessing intent.

Average scores across runs within each task, then average task means using declared weights. Report each family's task count, planned and completed runs, score, and interval, alongside the aggregate. Retain attempted failures under the declared scoring policy and report unattempted tasks separately; never silently drop either to improve a score.

For matched model comparisons, compute task-level paired score differences. Bootstrap task IDs with replacement, preserving all compared models for each sampled task; use a hierarchical within-task run resample when estimating run variability as well. Record the confidence level, bootstrap count, random seed, weighting, missing-run handling, and resampling scheme. Cluster related tasks by shared structure when dependence would otherwise exaggerate precision.

Distinguish uncertainty across tasks from variability across repeated runs. With a very small pilot, emphasize coverage and observed failures; do not present narrow-looking intervals as population evidence. With one run per task, state that within-task stochastic variability is unmeasured.

Report paired differences and intervals, not just rank order or overlap of separate confidence intervals. An interval containing zero leaves the difference unresolved; it does not establish equivalence. Label multiple pairwise comparisons as exploratory unless multiplicity was handled in advance. If ranking matters, report bootstrap leader frequencies and the treatment of ties.

Do not pool incompatible task subsets or attribute differences across separate grading passes solely to access changes. The source finding does not establish a universal ordering of model, harness, access, and task-composition effects.

#### 7. Validate the report

Before delivery, confirm that each scored numeric answer links to a file checksum, explicit selection, reproducible calculation, units, tolerance, and trial record. Reconcile all task and run counts against the manifest. State exclusions, incomplete work, unresolved selection ambiguity, and deviations from the planned protocol.

Completion means verified assets, at least one independently checked coordinate calculation, traceable trial results, and a report with explicit scope and limitations. If these conditions are unmet, label the output partial or blocked rather than successful.

## Inputs / Outputs

**Inputs**

- Task prompts and identifiers, with family labels or enough information to assign them.
- Protein coordinate files or retrievable structure references, plus intended chains, residues, models, and assembly.
- Verified reference answers or a specification for independent coordinate calculations.
- Agent endpoints or existing rollout logs; harness, tool, access, and budget settings.
- Scoring rules, tolerances, run count, and the comparison or diagnostic question.

**Outputs**

- An asset and selection manifest containing checksums, versions, task families, and structural conventions.
- A reference-calculation record containing selected atoms, methods, units, tolerances, and validation evidence.
- A trial ledger with one record per attempt: task ID, family, configuration ID, run ID, status, answer, score, error category, timing, and trace location.
- A comparison report with family-level and aggregate results, denominators, paired differences, uncertainty method, access checks, and limitations.
- A completion status of successful, partial, or blocked, with remaining work stated explicitly.

These are logical deliverables; return them inline or save them only within the output scope authorized for the evaluation.

## References

- [ProteinGeoBench: Frontier Structural-Biology Agents Are Fluent at Lookup but Lost in Geometry, arXiv:2610.26829v1](http://arxiv.org/abs/2610.26829v1). Research evidence for task-family separation, repeated-run uncertainty, and configuration-sensitive comparisons.
- [ProteinGeoBench repository](https://github.com/Protege-AI/ProteinGeoBench). Benchmark assets linked by the paper; verify the selected revision, task files, runner, and verifiers before execution.
