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

from agent_framework import ChatMessage, FunctionCallContent, FunctionResultContent, TextContent

from agent_framework_ag_ui._orchestration._message_hygiene import (
    deduplicate_messages,
    sanitize_tool_history,
)


def test_sanitize_tool_history_injects_confirm_changes_result() -> None:
    messages = [
        ChatMessage(
            role="assistant",
            contents=[
                FunctionCallContent(
                    name="confirm_changes",
                    call_id="call_confirm_123",
                    arguments='{"changes": "test"}',
                )
            ],
        ),
        ChatMessage(
            role="user",
            contents=[TextContent(text='{"accepted": true}')],
        ),
    ]

    sanitized = sanitize_tool_history(messages)

    tool_messages = [
        msg for msg in sanitized if (msg.role.value if hasattr(msg.role, "value") else str(msg.role)) == "tool"
    ]
    assert len(tool_messages) == 1
    assert str(tool_messages[0].contents[0].call_id) == "call_confirm_123"
    assert tool_messages[0].contents[0].result == "Confirmed"


def test_deduplicate_messages_prefers_non_empty_tool_results() -> None:
    messages = [
        ChatMessage(
            role="tool",
            contents=[FunctionResultContent(call_id="call1", result="")],
        ),
        ChatMessage(
            role="tool",
            contents=[FunctionResultContent(call_id="call1", result="result data")],
        ),
    ]

    deduped = deduplicate_messages(messages)
    assert len(deduped) == 1
    assert deduped[0].contents[0].result == "result data"

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
