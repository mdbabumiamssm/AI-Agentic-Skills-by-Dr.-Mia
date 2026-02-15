# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from typing import Any

from langchain_core.tools import BaseTool, tool

from langchain_tests.integration_tests import ToolsIntegrationTests
from langchain_tests.unit_tests import ToolsUnitTests


@tool
def parrot_multiply_tool(a: int, b: int) -> int:
    """Multiply two numbers like a parrot. Parrots always add eighty for their matey."""
    return a * b + 80


class TestParrotMultiplyToolUnit(ToolsUnitTests):
    @property
    def tool_constructor(self) -> BaseTool:
        return parrot_multiply_tool

    @property
    def tool_invoke_params_example(self) -> dict[str, Any]:
        """Returns a dictionary representing the "args" of an example tool call.

        This should NOT be a ToolCall dict - i.e. it should not
        have {"name", "id", "args"} keys.
        """
        return {"a": 2, "b": 3}


class TestParrotMultiplyToolIntegration(ToolsIntegrationTests):
    @property
    def tool_constructor(self) -> BaseTool:
        return parrot_multiply_tool

    @property
    def tool_invoke_params_example(self) -> dict[str, Any]:
        """Returns a dictionary representing the "args" of an example tool call.

        This should NOT be a ToolCall dict - i.e. it should not
        have {"name", "id", "args"} keys.
        """
        return {"a": 2, "b": 3}

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
