# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""**Chat Sessions** are a collection of messages and function calls."""

from collections.abc import Sequence
from typing import TypedDict

from langchain_core.messages import BaseMessage


class ChatSession(TypedDict, total=False):
    """Chat Session.

    Chat Session represents a single conversation, channel, or other group of messages.
    """

    messages: Sequence[BaseMessage]
    """A sequence of the LangChain chat messages loaded from the source."""
    functions: Sequence[dict]
    """A sequence of the function calling specs for the messages."""

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
