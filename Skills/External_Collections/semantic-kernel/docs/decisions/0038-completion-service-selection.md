<!--
# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

-->

---
# These are optional elements. Feel free to remove any of them.
status: accepted
contact: markwallace-microsoft
date: 2024-03-14
deciders: sergeymenshykh, markwallace, rbarreto, dmytrostruk
consulted: 
informed: 
---

# Completion Service Selection Strategy

## Context and Problem Statement

Today, SK uses the current `IAIServiceSelector` implementation to determine which type of service is used when running a text prompt.
The `IAIServiceSelector` implementation will return either a chat completion service, text generation service or it could return a service that implements both.
The prompt will be run using chat completion by default and falls back to text generation as the alternate option.

The behavior supersedes that description in [ADR-0015](0015-completion-service-selection.md)

## Decision Drivers

- Chat completion services are becoming dominant in the industry e.g. OpenAI has deprecated most of it's text generation services.
- Chat completion generally provides better responses and the ability to use advanced features e.g. tool calling.

## Decision Outcome

Chosen option: Keep the current behavior as described above.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->