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

import asyncio
from typing import Annotated

from semantic_kernel.exceptions import FunctionExecutionException
from semantic_kernel.functions.kernel_function_decorator import kernel_function
from semantic_kernel.kernel_pydantic import KernelBaseModel


class WaitPlugin(KernelBaseModel):
    """WaitPlugin provides a set of functions to wait for a certain amount of time.

    Usage:
        kernel.add_plugin(WaitPlugin(), plugin_name="wait")

    Examples:
        {{wait.wait 5}} => Wait for 5 seconds
    """

    @kernel_function
    async def wait(self, input: Annotated[float | str, "The number of seconds to wait, can be str or float."]) -> None:
        """Wait for a certain number of seconds."""
        if isinstance(input, str):
            try:
                input = float(input)
            except ValueError as exc:
                raise FunctionExecutionException("seconds text must be a number") from exc
        await asyncio.sleep(abs(input))

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
