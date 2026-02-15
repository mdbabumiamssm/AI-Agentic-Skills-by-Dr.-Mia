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

# pathway-analysis

## Overview

Functional enrichment and pathway analysis using R/Bioconductor. Supports over-representation analysis (ORA) and Gene Set Enrichment Analysis (GSEA) across Gene Ontology, KEGG, Reactome, and WikiPathways databases.

**Tool type:** r | **Primary tools:** clusterProfiler, ReactomePA, rWikiPathways, enrichplot

## Skills

| Skill | Description |
|-------|-------------|
| go-enrichment | Gene Ontology over-representation analysis with enrichGO |
| kegg-pathways | KEGG pathway enrichment with enrichKEGG and enrichMKEGG |
| reactome-pathways | Reactome pathway enrichment with ReactomePA |
| wikipathways | WikiPathways enrichment with enrichWP and rWikiPathways |
| gsea | Gene Set Enrichment Analysis with gseGO, gseKEGG |
| enrichment-visualization | Dot plots, bar plots, enrichment maps, cnetplots, GSEA plots |

## Example Prompts

- "Run GO enrichment on my differentially expressed genes"
- "Find enriched biological processes for these genes"
- "What molecular functions are over-represented in my gene list?"
- "Find enriched KEGG pathways for my gene set"
- "What pathways are active in my differentially expressed genes?"
- "Run KEGG module enrichment analysis"
- "Run Reactome pathway enrichment on my genes"
- "Find enriched Reactome pathways for my DEGs"
- "Run WikiPathways enrichment analysis"
- "Find community-curated pathways for my gene list"
- "Run GSEA on my ranked gene list"
- "Perform gene set enrichment analysis using GO terms"
- "Run GSEA with KEGG pathways"
- "Create a dot plot of my enrichment results"
- "Make an enrichment map showing term relationships"
- "Show a gene-concept network for top pathways"
- "Create a GSEA running score plot"

## Requirements

```r
BiocManager::install(c('clusterProfiler', 'enrichplot', 'org.Hs.eg.db'))
BiocManager::install(c('ReactomePA', 'rWikiPathways'))
```

## Related Skills

- **differential-expression** - Generate gene lists and statistics for enrichment
- **single-cell** - Marker genes can be analyzed with pathway enrichment
- **database-access** - Fetch gene annotations from NCBI


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->