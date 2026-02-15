# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

# Copyright (c) Microsoft. All rights reserved.

"""GAIA benchmark module for Agent Framework."""

import importlib.metadata

from ._types import Evaluation, Evaluator, Prediction, Task, TaskResult, TaskRunner
from .gaia import GAIA, GAIATelemetryConfig, gaia_scorer, viewer_main

try:
    __version__ = importlib.metadata.version(__name__)
except importlib.metadata.PackageNotFoundError:
    __version__ = "0.0.0"  # Fallback for development mode

__all__ = [
    "GAIA",
    "Evaluation",
    "Evaluator",
    "GAIATelemetryConfig",
    "Prediction",
    "Task",
    "TaskResult",
    "TaskRunner",
    "gaia_scorer",
    "viewer_main",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
