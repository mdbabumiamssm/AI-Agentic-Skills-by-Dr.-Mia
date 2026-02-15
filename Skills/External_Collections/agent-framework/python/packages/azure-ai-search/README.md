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

# Get Started with Microsoft Agent Framework Azure AI Search

Please install this package via pip:

```bash
pip install agent-framework-azure-ai-search --pre
```

## Azure AI Search Integration

The Azure AI Search integration provides context providers for RAG (Retrieval Augmented Generation) capabilities with two modes:

- **Semantic Mode**: Fast hybrid search (vector + keyword) with semantic ranking
- **Agentic Mode**: Multi-hop reasoning using Knowledge Bases for complex queries

### Basic Usage Example

See the [Azure AI Search context provider examples](https://github.com/microsoft/agent-framework/tree/main/python/samples/getting_started/agents/azure_ai/) which demonstrate:

- Semantic search with hybrid (vector + keyword) queries
- Agentic mode with Knowledge Bases for complex multi-hop reasoning
- Environment variable configuration with Settings class
- API key and managed identity authentication


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->