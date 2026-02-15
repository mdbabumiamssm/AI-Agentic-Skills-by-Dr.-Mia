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

# Pharmacovigilance Agent

**ID:** `biomedical.clinical.pharmacovigilance`
**Version:** 1.0.0
**Status:** Production
**Category:** Clinical / Safety

---

## Overview

The **Pharmacovigilance Agent** monitors real-world data sources to detect potential Adverse Drug Events (ADEs) that were not caught during clinical trials. It combines structured database analysis (FDA FAERS) with unstructured social media mining.

## Key Capabilities

### 1. Social Media Mining
- **NER:** Uses BERT-based models to extract drug names and symptoms from Twitter/Reddit posts.
- **Sentiment Analysis:** Filters for negative experiences.
- **Signal Detection:** "Is there a spike in 'Headache' mentions associated with 'Drug X' this week?"

### 2. FAERS Analysis
- **Data Source:** FDA Adverse Event Reporting System.
- **Statistics:** Calculates Proportional Reporting Ratios (PRR) and Reporting Odds Ratios (ROR) to identify statistical signals.

## Ethics
*Strictly adheres to data privacy. Social media data is anonymized and aggregated.*

## References
- *Web-RADR*
- *SIDER Side Effect Resource*


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->