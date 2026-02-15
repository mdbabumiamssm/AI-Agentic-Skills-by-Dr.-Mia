# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""__ModuleName__ tools."""

from typing import Type

from langchain_core.callbacks import (
    CallbackManagerForToolRun,
)
from langchain_core.tools import BaseTool
from pydantic import BaseModel, Field


class __ModuleName__ToolInput(BaseModel):
    """Input schema for __ModuleName__ tool.

    This docstring is **not** part of what is sent to the model when performing tool
    calling. The Field default values and descriptions **are** part of what is sent to
    the model when performing tool calling.
    """

    # TODO: Add input args and descriptions.
    a: int = Field(..., description="first number to add")
    b: int = Field(..., description="second number to add")


class __ModuleName__Tool(BaseTool):  # type: ignore[override]
    """__ModuleName__ tool.

    Setup:
        # TODO: Replace with relevant packages, env vars.
        Install `__package_name__` and set environment variable
        `__MODULE_NAME___API_KEY`.

        ```bash
        pip install -U __package_name__
        export __MODULE_NAME___API_KEY="your-api-key"
        ```

    Instantiation:
        ```python
        tool = __ModuleName__Tool(
            # TODO: init params
        )
        ```

    Invocation with args:
        ```python
        # TODO: invoke args
        tool.invoke({...})
        ```

        ```python
        # TODO: output of invocation
        ```

    Invocation with ToolCall:

        ```python
        # TODO: invoke args
        tool.invoke({"args": {...}, "id": "1", "name": tool.name, "type": "tool_call"})
        ```

        ```python
        # TODO: output of invocation

        ```
    """  # noqa: E501

    # TODO: Set tool name and description
    name: str = "TODO: Tool name"
    """The name that is passed to the model when performing tool calling."""
    description: str = "TODO: Tool description."
    """The description that is passed to the model when performing tool calling."""
    args_schema: Type[BaseModel] = __ModuleName__ToolInput
    """The schema that is passed to the model when performing tool calling."""

    # TODO: Add any other init params for the tool.
    # param1: str | None
    # """param1 determines foobar"""

    # TODO: Replaced (a, b) with real tool arguments.
    def _run(
        self, a: int, b: int, *, run_manager: CallbackManagerForToolRun | None = None
    ) -> str:
        return str(a + b + 80)

    # TODO: Implement if tool has native async functionality, otherwise delete.

    # async def _arun(
    #     self,
    #     a: int,
    #     b: int,
    #     *,
    #     run_manager: AsyncCallbackManagerForToolRun | None = None,
    # ) -> str:
    #     ...

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
