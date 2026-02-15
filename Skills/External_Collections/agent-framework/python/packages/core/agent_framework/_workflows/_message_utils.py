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

"""Shared helpers for normalizing workflow message inputs."""

from collections.abc import Sequence

from agent_framework import ChatMessage, Role


def normalize_messages_input(
    messages: str | ChatMessage | Sequence[str | ChatMessage] | None = None,
) -> list[ChatMessage]:
    """Normalize heterogeneous message inputs to a list of ChatMessage objects.

    Args:
        messages: String, ChatMessage, or sequence of either. None yields empty list.

    Returns:
        List of ChatMessage instances suitable for workflow consumption.
    """
    if messages is None:
        return []

    if isinstance(messages, str):
        return [ChatMessage(role=Role.USER, text=messages)]

    if isinstance(messages, ChatMessage):
        return [messages]

    normalized: list[ChatMessage] = []
    for item in messages:
        if isinstance(item, str):
            normalized.append(ChatMessage(role=Role.USER, text=item))
        elif isinstance(item, ChatMessage):
            normalized.append(item)
        else:
            raise TypeError(
                f"Messages sequence must contain only str or ChatMessage instances; found {type(item).__name__}."
            )
    return normalized


__all__ = ["normalize_messages_input"]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
