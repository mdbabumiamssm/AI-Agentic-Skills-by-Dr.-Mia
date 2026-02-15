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

# PyLabRobot Control Skill

**ID:** `biomedical.lab_automation.pylabrobot`
**Version:** 1.0.0
**Status:** Alpha
**Category:** Lab Automation / Hardware Agnostic

---

## Overview

**PyLabRobot** is a hardware-agnostic Python library that allows a single agent to control multiple types of liquid handling robots (Hamilton, Tecan, Opentrons) using a unified API. This skill enables the "Universal Lab Agent" concept.

---

## Capabilities

- **Unified Commands:** Use `await lh.aspirate()` regardless of the underlying hardware.
- **Visualizer:** Web-based simulation of deck states before execution.
- **Resource Management:** Tracks volume and liquid classes (viscosity, volatility).

## Integration

Use this skill when the specific robot model is unknown or when building a facility-wide orchestration agent.

## References
- [PyLabRobot GitHub](https://github.com/PyLabRobot/pylabrobot)


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->