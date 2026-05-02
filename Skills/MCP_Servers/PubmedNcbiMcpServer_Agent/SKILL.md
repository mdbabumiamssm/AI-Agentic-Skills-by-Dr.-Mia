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
name: 'pubmed-ncbi-mcp-server'
description: 'MCP server exposing NCBI E-utilities for PubMed search, article metadata and full-text retrieval, citation generation, MeSH exploration, and related-article discovery via STDIO or Streamable HTTP.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# PubMed NCBI E-utilities MCP Server

## Overview
This skill wires an AI agent into the cyanheads `pubmed-mcp-server`, a TypeScript MCP server that wraps NCBI's E-utilities (ESearch, ESummary, EFetch, ELink) along with MeSH and Unpaywall lookups. It enables agents to perform reproducible literature mining — searching PubMed, fetching metadata and open-access full text, formatting citations, and discovering related research — over either STDIO or Streamable HTTP transports.

## When to Use This Skill
- Building a clinical or biomedical research agent that must ground answers in primary literature.
- Performing PubMed searches with controlled MeSH vocabulary or PMID/DOI lookups inside an agent loop.
- Retrieving article metadata, abstracts, or open-access full text for downstream summarization, RAG, or evidence ranking.
- Generating formatted citations (e.g., for systematic reviews, scoping reviews, or report drafting).
- Discovering related articles via NCBI ELink to expand a literature corpus.
- Hosting a shared MCP endpoint (Streamable HTTP) that multiple agents or LLM clients can query.

## Core Capabilities
1. **PubMed search (ESearch)** — Run keyword or fielded queries against PubMed and return PMIDs with optional pagination, sort, and date filters.
2. **Article metadata (ESummary / EFetch)** — Resolve PMIDs to structured metadata: title, authors, journal, year, abstract, MeSH headings, and identifiers (DOI, PMC).
3. **Full-text retrieval** — Fetch open-access full text where available via PMC and Unpaywall, returning XML/text suitable for downstream parsing.
4. **Citation generation** — Format article records into common citation styles for inclusion in agent outputs and reports.
5. **MeSH exploration** — Look up MeSH terms and tree numbers to refine search strategies with controlled vocabulary.
6. **Related-article discovery (ELink)** — Surface PubMed-related neighbors of a given PMID to expand or seed a corpus.
7. **Article-ID conversion** — Translate between PMID, PMCID, and DOI to stitch records across NCBI and publisher systems.
8. **Dual transport** — Run as a local STDIO MCP server for desktop clients (Claude Desktop, IDEs) or as a Streamable HTTP server for shared/remote agent access.

## Inputs / Outputs
**Inputs**
- Natural-language or fielded PubMed queries (e.g., `"GLP-1[Title] AND 2024[PDAT]"`).
- Article identifiers: PMIDs, PMCIDs, or DOIs.
- MeSH terms or tree-number lookups.
- Optional NCBI API key (recommended) and Unpaywall email for open-access resolution; transport selection (STDIO or HTTP) and port configuration.

**Outputs**
- Lists of PMIDs and pagination cursors from search.
- Structured article metadata records (JSON) with abstracts and MeSH terms.
- Full-text payloads (XML/plain text) when open access is available.
- Formatted citation strings.
- Related-article PMID sets and ID-conversion maps.

## References
- Source repository: https://github.com/cyanheads/pubmed-mcp-server
- NCBI E-utilities documentation: https://www.ncbi.nlm.nih.gov/books/NBK25501/
- PubMed: https://pubmed.ncbi.nlm.nih.gov/
- PMC (PubMed Central) Open Access Subset: https://www.ncbi.nlm.nih.gov/pmc/tools/openftlist/
- MeSH (Medical Subject Headings): https://www.nlm.nih.gov/mesh/meshhome.html
- Unpaywall API: https://unpaywall.org/products/api
- Model Context Protocol specification: https://modelcontextprotocol.io
