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

import pytest

from semantic_kernel.connectors.ai.open_ai.services.open_ai_chat_completion import OpenAIChatCompletion
from semantic_kernel.kernel import Kernel


def test_kernel_deep_copy_fail_with_services():
    kernel = Kernel()
    kernel.add_service(OpenAIChatCompletion())

    with pytest.raises(TypeError):
        # This will fail because OpenAIChatCompletion is not serializable, more specifically,
        # the client is not serializable
        kernel.model_copy(deep=True)


async def test_kernel_clone():
    kernel = Kernel()
    kernel.add_service(OpenAIChatCompletion())

    kernel_clone = kernel.clone()

    assert kernel_clone is not None
    assert kernel_clone.services is not None and len(kernel_clone.services) > 0

    function_result = await kernel.invoke_prompt("Hello World")
    assert function_result is not None
    assert function_result.value is not None
    assert len(str(function_result)) > 0

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
