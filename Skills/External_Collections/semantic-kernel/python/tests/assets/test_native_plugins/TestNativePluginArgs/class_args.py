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

from typing import Annotated

from semantic_kernel.functions.kernel_function_decorator import kernel_function


class TestNativeEchoBotPlugin:
    """Description: Test Native Plugin for testing purposes"""

    def __init__(self, static_input: str | None = None):
        self.static_input = static_input or ""

    @kernel_function(
        description="Echo for input text with static",
        name="echo",
    )
    def echo(self, text: Annotated[str, "The text to echo"]) -> str:
        """Echo for input text with a static input

        Example:
            "hello world" => "hello world"
        Args:
            text -- The text to echo

        Returns:
            input text
        """
        return self.static_input + text

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
