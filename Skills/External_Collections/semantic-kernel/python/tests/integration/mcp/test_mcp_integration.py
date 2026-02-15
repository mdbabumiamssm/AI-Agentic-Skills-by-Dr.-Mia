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


import os
from typing import TYPE_CHECKING

from semantic_kernel.connectors.mcp import MCPStdioPlugin
from semantic_kernel.functions.kernel_arguments import KernelArguments

if TYPE_CHECKING:
    from semantic_kernel import Kernel


async def test_from_mcp(kernel: "Kernel"):
    mcp_server_path = os.path.join(
        os.path.dirname(__file__), "../../assets/test_plugins", "TestMCPPlugin", "mcp_server.py"
    )
    async with MCPStdioPlugin(
        name="TestMCPPlugin",
        command="python",
        args=[mcp_server_path],
    ) as plugin:
        assert plugin is not None
        assert plugin.name == "TestMCPPlugin"

        loaded_plugin = kernel.add_plugin(plugin)

        assert loaded_plugin is not None
        assert loaded_plugin.name == "TestMCPPlugin"
        assert len(loaded_plugin.functions) == 2

        result = await loaded_plugin.functions["echo_tool"].invoke(kernel, arguments=KernelArguments(message="test"))
        assert "Tool echo: test" in result.value[0].text

        result = await loaded_plugin.functions["echo_prompt"].invoke(kernel, arguments=KernelArguments(message="test"))
        assert "test" in result.value[0].content

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
