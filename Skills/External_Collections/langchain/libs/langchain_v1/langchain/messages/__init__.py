# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Message and message content types.

Includes message types for different roles (e.g., human, AI, system), as well as types
for message content blocks (e.g., text, image, audio) and tool calls.
"""

from langchain_core.messages import (
    AIMessage,
    AIMessageChunk,
    Annotation,
    AnyMessage,
    AudioContentBlock,
    Citation,
    ContentBlock,
    DataContentBlock,
    FileContentBlock,
    HumanMessage,
    ImageContentBlock,
    InputTokenDetails,
    InvalidToolCall,
    MessageLikeRepresentation,
    NonStandardAnnotation,
    NonStandardContentBlock,
    OutputTokenDetails,
    PlainTextContentBlock,
    ReasoningContentBlock,
    RemoveMessage,
    ServerToolCall,
    ServerToolCallChunk,
    ServerToolResult,
    SystemMessage,
    TextContentBlock,
    ToolCall,
    ToolCallChunk,
    ToolMessage,
    UsageMetadata,
    VideoContentBlock,
    trim_messages,
)

__all__ = [
    "AIMessage",
    "AIMessageChunk",
    "Annotation",
    "AnyMessage",
    "AudioContentBlock",
    "Citation",
    "ContentBlock",
    "DataContentBlock",
    "FileContentBlock",
    "HumanMessage",
    "ImageContentBlock",
    "InputTokenDetails",
    "InvalidToolCall",
    "MessageLikeRepresentation",
    "NonStandardAnnotation",
    "NonStandardContentBlock",
    "OutputTokenDetails",
    "PlainTextContentBlock",
    "ReasoningContentBlock",
    "RemoveMessage",
    "ServerToolCall",
    "ServerToolCallChunk",
    "ServerToolResult",
    "SystemMessage",
    "TextContentBlock",
    "ToolCall",
    "ToolCallChunk",
    "ToolMessage",
    "UsageMetadata",
    "VideoContentBlock",
    "trim_messages",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
