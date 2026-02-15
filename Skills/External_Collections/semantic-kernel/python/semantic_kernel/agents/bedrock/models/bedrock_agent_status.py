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


from enum import Enum

from semantic_kernel.utils.feature_stage_decorator import experimental


@experimental
class BedrockAgentStatus(str, Enum):
    """Bedrock Agent Status.

    https://docs.aws.amazon.com/bedrock/latest/APIReference/API_agent_PrepareAgent.html#API_agent_PrepareAgent_ResponseElements
    """

    CREATING = "CREATING"
    PREPARING = "PREPARING"
    PREPARED = "PREPARED"
    NOT_PREPARED = "NOT_PREPARED"
    DELETING = "DELETING"
    FAILED = "FAILED"
    VERSIONING = "VERSIONING"
    UPDATING = "UPDATING"

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
