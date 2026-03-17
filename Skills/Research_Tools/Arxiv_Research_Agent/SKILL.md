# Arxiv Research Agent Skill

**Domain:** Research Tools / Academic Literature
**Status:** Active

## Overview
The Arxiv Research Agent allows AI agents to autonomously query the arXiv API to find the latest preprints and research papers. This skill is crucial for "Deep Research" workflows where agents need to compile state-of-the-art literature.

## Capabilities
- `search_papers(query, max_results)`: Searches the arXiv database using a keyword query and returns parsed metadata including titles, authors, links, and summaries.

## Implementation Details
- **Language:** Python
- **Dependencies:** Standard library (`urllib`, `xml.etree.ElementTree`)
- **Main Class:** `ArxivResearchAgent`
