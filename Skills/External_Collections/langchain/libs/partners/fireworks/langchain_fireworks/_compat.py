# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Converts between AIMessage output formats, governed by `output_version`."""

from __future__ import annotations

from langchain_core.messages import AIMessage


def _convert_from_v1_to_chat_completions(message: AIMessage) -> AIMessage:
    """Convert a v1 message to the Chat Completions format."""
    if isinstance(message.content, list):
        new_content: list = []
        for block in message.content:
            if isinstance(block, dict):
                block_type = block.get("type")
                if block_type == "text":
                    # Strip annotations
                    new_content.append({"type": "text", "text": block["text"]})
                elif block_type in ("reasoning", "tool_call"):
                    pass
                else:
                    new_content.append(block)
            else:
                new_content.append(block)
        return message.model_copy(update={"content": new_content})

    return message

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
