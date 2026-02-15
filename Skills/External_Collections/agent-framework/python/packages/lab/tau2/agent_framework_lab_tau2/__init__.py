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

"""Tau2 Benchmark for Agent Framework."""

import importlib.metadata

from ._tau2_utils import patch_env_set_state, unpatch_env_set_state
from .runner import ASSISTANT_AGENT_ID, ORCHESTRATOR_ID, USER_SIMULATOR_ID, TaskRunner

try:
    __version__ = importlib.metadata.version(__name__)
except importlib.metadata.PackageNotFoundError:
    __version__ = "0.0.0"  # Fallback for development mode

__all__ = [
    "ASSISTANT_AGENT_ID",
    "ORCHESTRATOR_ID",
    "USER_SIMULATOR_ID",
    "TaskRunner",
    "patch_env_set_state",
    "unpatch_env_set_state",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
