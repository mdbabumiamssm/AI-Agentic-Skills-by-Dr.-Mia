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

import logging
from dataclasses import dataclass

from semantic_kernel.contents.streaming_chat_message_content import StreamingChatMessageContent
from semantic_kernel.utils.feature_stage_decorator import experimental

logger: logging.Logger = logging.getLogger(__name__)


@experimental
@dataclass
class FunctionActionResult:
    """Function Action Result."""

    function_call_streaming_content: StreamingChatMessageContent
    function_result_streaming_content: StreamingChatMessageContent
    tool_outputs: list[dict[str, str]]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
