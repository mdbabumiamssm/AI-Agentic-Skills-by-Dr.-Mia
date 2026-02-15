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

# Medaka Polishing - Usage Guide

## Overview
Medaka uses neural networks to polish consensus sequences and call variants from Oxford Nanopore data. Models are trained on specific basecaller versions for optimal accuracy.

## Prerequisites
```bash
conda create -n medaka -c conda-forge -c bioconda medaka
conda activate medaka
```

## Quick Start
Tell your AI agent what you want to do:
- "Polish my draft assembly with medaka"
- "Call variants from ONT reads using medaka"

## Example Prompts

### Assembly Polishing
> "Polish my draft assembly draft.fa using reads.fastq.gz with medaka"

> "Run medaka consensus on my Nanopore assembly using the R10.4.1 SUP model"

### Variant Calling
> "Call haploid variants from my ONT reads aligned to the reference"

### Model Selection
> "What medaka model should I use for R10.4.1 flowcell with SUP basecalling?"

> "List available medaka models and help me choose the right one"

## What the Agent Will Do
1. Identify the correct medaka model based on your basecaller and chemistry
2. Run medaka_consensus for polishing or medaka_variant for variant calling
3. Handle intermediate files and cleanup
4. Report improvement in consensus quality

## Model Selection
Critical: Use the model matching your basecaller.

```bash
# List available models
medaka tools list_models

# Common choices:
# - r1041_e82_400bps_sup_v5.0.0 (R10.4.1, SUP, latest)
# - r941_min_sup_g507 (R9.4.1, SUP)
```

## Expected Accuracy

| Basecaller | Before Polish | After Polish |
|------------|---------------|--------------|
| SUP R10.4.1 | Q20 (~99%) | Q40+ (>99.99%) |
| HAC R10.4.1 | Q15-18 | Q30+ |

## Tips
- Wrong model is the most common cause of poor polishing - always check basecaller version
- For diploid variant calling, use Clair3 instead (medaka diploid is deprecated)
- Reduce batch size (`-b 50`) if running out of memory
- Process by region (`--region chr1`) for large genomes
- Use GPU if available for much faster processing

## Resources
- [Medaka GitHub](https://github.com/nanoporetech/medaka)
- [ONT Community](https://community.nanoporetech.com/)


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->