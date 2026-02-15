# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""__ModuleName__ toolkits."""

from typing import List

from langchain_core.tools import BaseTool, BaseToolkit


class __ModuleName__Toolkit(BaseToolkit):
    # TODO: Replace all TODOs in docstring. See example docstring:
    # https://github.com/langchain-ai/langchain/blob/c123cb2b304f52ab65db4714eeec46af69a861ec/libs/community/langchain_community/agent_toolkits/sql/toolkit.py#L19
    """__ModuleName__ toolkit.

    # TODO: Replace with relevant packages, env vars, etc.
    Setup:
        Install `__package_name__` and set environment variable
        `__MODULE_NAME___API_KEY`.

        ```bash
        pip install -U __package_name__
        export __MODULE_NAME___API_KEY="your-api-key"
        ```

    # TODO: Populate with relevant params.
    Key init args:
        arg 1: type
            description
        arg 2: type
            description

    # TODO: Replace with relevant init params.
    Instantiate:
        ```python
        from __package_name__ import __ModuleName__Toolkit

        toolkit = __ModuleName__Toolkit(
            # ...
        )
        ```

    Tools:
        ```python
        toolkit.get_tools()
        ```

        ```txt
        # TODO: Example output.
        ```

    Use within an agent:
        ```python
        from langgraph.prebuilt import create_react_agent

        agent_executor = create_react_agent(llm, tools)

        example_query = "..."

        events = agent_executor.stream(
            {"messages": [("user", example_query)]},
            stream_mode="values",
        )
        for event in events:
            event["messages"][-1].pretty_print()
        ```

        ```txt
        # TODO: Example output.
        ```

    """

    # TODO: This method must be implemented to list tools.
    def get_tools(self) -> List[BaseTool]:
        raise NotImplementedError()

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
