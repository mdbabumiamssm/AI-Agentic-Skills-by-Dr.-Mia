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

# Example: Executive Review Prep

**User Request**: "Prep for quarterly executive review on Friday"

## Workflow

### 1. Search for Context
```
Notion:notion-search
query: "Q4 objectives" + "KPIs" + "quarterly results"
```

Found:
- Q4 OKRs and progress
- Product metrics dashboard
- Engineering velocity reports
- Customer feedback summary

### 2. Fetch & Analyze
```
Notion:notion-fetch (5 pages)
```

**Key metrics**:
- **Revenue**: $2.4M ARR (96% of Q4 target)
- **Customer Growth**: 145 new customers (exceeds 120 target)
- **Churn**: 3.2% (below 5% target)
- **Product**: 3 major features shipped, 2 in beta
- **Engineering**: 94% uptime (above 95% SLA)

### 3. Add Codex Research Context
Added context on:
- Industry benchmarks for SaaS metrics
- Typical Q4 sales patterns
- Best practices for executive presentations

### 4. Create Pre-Read (Internal)
```
Notion:notion-create-pages
title: "Q4 Review - Pre-Read (Internal)"
```

**Pre-read sections**:
- **Executive Summary**: Strong quarter, missed revenue by 4% but exceeded customer growth
- **Detailed Metrics**: All KPIs with trend lines
- **Wins**: Product launches, key customer acquisitions
- **Challenges**: Sales pipeline conversion, engineering hiring
- **Q1 Preview**: Strategic priorities

### 5. Create Presentation Agenda
```
Notion:notion-create-pages
title: "Q4 Executive Review - Agenda"
```

**Agenda** (90 min):
- Q4 Results Overview (15 min)
- Revenue & Growth Deep Dive (20 min)
- Product & Engineering Update (20 min)
- Customer Success Highlights (15 min)
- Q1 Strategic Plan (15 min)
- Discussion & Questions (15 min)

### 6. Link Supporting Docs
Connected to OKRs, metrics dashboards, and Q1 planning docs.

## Outputs

**Internal Pre-Read**: Comprehensive context with honest assessment
**Executive Agenda**: Structured 90-min presentation
**Both in Notion** with links to supporting data

## Key Success Factors
- Synthesized data from multiple sources (OKRs, metrics, feedback)
- Added industry context and benchmarks
- Created honest internal assessment (not just wins)
- Structured agenda with time allocations
- Linked to source data for drill-down during Q&A


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->