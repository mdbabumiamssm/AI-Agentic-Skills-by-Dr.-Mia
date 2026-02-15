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

from dataclasses import dataclass

from semantic_kernel.agents.runtime.core.agent_id import AgentId
from semantic_kernel.agents.runtime.core.cancellation_token import CancellationToken
from semantic_kernel.agents.runtime.core.topic import TopicId
from semantic_kernel.utils.feature_stage_decorator import experimental


@experimental
@dataclass
class MessageContext:
    """Context for a message sent to an agent."""

    sender: AgentId | None
    topic_id: TopicId | None
    is_rpc: bool
    cancellation_token: CancellationToken
    message_id: str

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
