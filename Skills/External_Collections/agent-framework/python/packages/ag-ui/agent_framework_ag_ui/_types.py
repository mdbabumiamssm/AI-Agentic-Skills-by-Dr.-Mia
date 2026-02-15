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

"""Type definitions for AG-UI integration."""

from typing import Any, TypedDict


class PredictStateConfig(TypedDict):
    """Configuration for predictive state updates."""

    state_key: str
    tool: str
    tool_argument: str | None


class RunMetadata(TypedDict):
    """Metadata for agent run."""

    run_id: str
    thread_id: str
    predict_state: list[PredictStateConfig] | None


class AgentState(TypedDict):
    """Base state for AG-UI agents."""

    messages: list[Any] | None

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
