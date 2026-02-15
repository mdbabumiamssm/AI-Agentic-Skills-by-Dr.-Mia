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

This page is a list of issues you could encounter along with their fixes.

# Forge
**Poetry configuration invalid**

The poetry configuration is invalid: 
- Additional properties are not allowed ('group' was unexpected)
<img width="487" alt="Screenshot 2023-09-22 at 5 42 59 PM" src="https://github.com/Significant-Gravitas/AutoGPT/assets/9652976/dd451e6b-8114-44de-9928-075f5f06d661">

**Pydantic Validation Error**

Remove your sqlite agent.db file. it's probably because some of your data is not complying with the new spec (we will create migrations soon to avoid this problem)


*Solution*

Update poetry

# Benchmark
TODO

# Frontend
TODO


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->