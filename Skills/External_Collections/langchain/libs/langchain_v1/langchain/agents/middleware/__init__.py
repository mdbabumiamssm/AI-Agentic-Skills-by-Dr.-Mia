# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Entrypoint to using [middleware](https://docs.langchain.com/oss/python/langchain/middleware) plugins with [Agents](https://docs.langchain.com/oss/python/langchain/agents)."""  # noqa: E501

from .context_editing import (
    ClearToolUsesEdit,
    ContextEditingMiddleware,
)
from .file_search import FilesystemFileSearchMiddleware
from .human_in_the_loop import (
    HumanInTheLoopMiddleware,
    InterruptOnConfig,
)
from .model_call_limit import ModelCallLimitMiddleware
from .model_fallback import ModelFallbackMiddleware
from .model_retry import ModelRetryMiddleware
from .pii import PIIDetectionError, PIIMiddleware
from .shell_tool import (
    CodexSandboxExecutionPolicy,
    DockerExecutionPolicy,
    HostExecutionPolicy,
    RedactionRule,
    ShellToolMiddleware,
)
from .summarization import SummarizationMiddleware
from .todo import TodoListMiddleware
from .tool_call_limit import ToolCallLimitMiddleware
from .tool_emulator import LLMToolEmulator
from .tool_retry import ToolRetryMiddleware
from .tool_selection import LLMToolSelectorMiddleware
from .types import (
    AgentMiddleware,
    AgentState,
    ModelRequest,
    ModelResponse,
    after_agent,
    after_model,
    before_agent,
    before_model,
    dynamic_prompt,
    hook_config,
    wrap_model_call,
    wrap_tool_call,
)

__all__ = [
    "AgentMiddleware",
    "AgentState",
    "ClearToolUsesEdit",
    "CodexSandboxExecutionPolicy",
    "ContextEditingMiddleware",
    "DockerExecutionPolicy",
    "FilesystemFileSearchMiddleware",
    "HostExecutionPolicy",
    "HumanInTheLoopMiddleware",
    "InterruptOnConfig",
    "LLMToolEmulator",
    "LLMToolSelectorMiddleware",
    "ModelCallLimitMiddleware",
    "ModelFallbackMiddleware",
    "ModelRequest",
    "ModelResponse",
    "ModelRetryMiddleware",
    "PIIDetectionError",
    "PIIMiddleware",
    "RedactionRule",
    "ShellToolMiddleware",
    "SummarizationMiddleware",
    "TodoListMiddleware",
    "ToolCallLimitMiddleware",
    "ToolRetryMiddleware",
    "after_agent",
    "after_model",
    "before_agent",
    "before_model",
    "dynamic_prompt",
    "hook_config",
    "wrap_model_call",
    "wrap_tool_call",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
