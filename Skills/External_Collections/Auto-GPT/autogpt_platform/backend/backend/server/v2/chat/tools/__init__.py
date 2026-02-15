# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from typing import TYPE_CHECKING, Any

from openai.types.chat import ChatCompletionToolParam

from backend.server.v2.chat.model import ChatSession

from .base import BaseTool
from .find_agent import FindAgentTool
from .run_agent import RunAgentTool

if TYPE_CHECKING:
    from backend.server.v2.chat.response_model import StreamToolExecutionResult

# Initialize tool instances
find_agent_tool = FindAgentTool()
run_agent_tool = RunAgentTool()

# Export tools as OpenAI format
tools: list[ChatCompletionToolParam] = [
    find_agent_tool.as_openai_tool(),
    run_agent_tool.as_openai_tool(),
]


async def execute_tool(
    tool_name: str,
    parameters: dict[str, Any],
    user_id: str | None,
    session: ChatSession,
    tool_call_id: str,
) -> "StreamToolExecutionResult":

    tool_map: dict[str, BaseTool] = {
        "find_agent": find_agent_tool,
        "run_agent": run_agent_tool,
    }
    if tool_name not in tool_map:
        raise ValueError(f"Tool {tool_name} not found")
    return await tool_map[tool_name].execute(
        user_id, session, tool_call_id, **parameters
    )

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
