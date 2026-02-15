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

# Scanpy Benchmarks

This directory contains code for benchmarking Scanpy using [asv][].

The functionality is checked using the [`benchmark.yml`][] workflow.
Benchmarks are run using the [benchmark bot][].

[asv]: https://asv.readthedocs.io/
[`benchmark.yml`]: ../.github/workflows/benchmark.yml
[benchmark bot]: https://github.com/apps/scverse-benchmark

## Data processing in benchmarks

Each dataset is processed so it has

- `.layers['counts']` (containing data in C/row-major format) and `.layers['counts-off-axis']` (containing data in FORTRAN/column-major format)
- `.X` and `.layers['off-axis']` with log-transformed data (formats like above)
- a `.var['mt']` boolean column indicating mitochondrial genes

The benchmarks are set up so the `layer` parameter indicates the layer that will be moved into `.X` before the benchmark.
That way, we don’t need to add `layer=layer` everywhere.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->