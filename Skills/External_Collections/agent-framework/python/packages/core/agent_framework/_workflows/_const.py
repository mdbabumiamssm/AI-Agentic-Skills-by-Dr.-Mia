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

# Default maximum iterations for workflow execution.
DEFAULT_MAX_ITERATIONS = 100

# Key used to store executor state in shared state.
EXECUTOR_STATE_KEY = "_executor_state"

# Source identifier for internal workflow messages.
INTERNAL_SOURCE_PREFIX = "internal"

# SharedState key for storing run kwargs that should be passed to agent invocations.
# Used by all orchestration patterns (Sequential, Concurrent, GroupChat, Handoff, Magentic)
# to pass kwargs from workflow.run_stream() through to agent.run_stream() and @ai_function tools.
WORKFLOW_RUN_KWARGS_KEY = "_workflow_run_kwargs"


def INTERNAL_SOURCE_ID(executor_id: str) -> str:
    """Generate an internal source ID for a given executor."""
    return f"{INTERNAL_SOURCE_PREFIX}:{executor_id}"

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
