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


from semantic_kernel import Kernel
from semantic_kernel.functions.function_result import FunctionResult
from semantic_kernel.functions.kernel_arguments import KernelArguments
from semantic_kernel.functions.kernel_function import KernelFunction
from semantic_kernel.functions.kernel_function_decorator import kernel_function
from semantic_kernel.text import aggregate_chunked_results


async def test_aggregate_results():
    kernel = Kernel()

    @kernel_function(name="func")
    def function(kernel, arguments):
        return FunctionResult(
            function=func.metadata,
            value=arguments["input"],
            metadata={},
        )

    func = KernelFunction.from_method(method=function, plugin_name="test")

    chunked = [
        "This is a test of the emergency broadcast system.",
        "This is only a test",
        "We repeat, this is only a test? A unit test",
        "A small note! And another? And once again!",
        "Seriously, this is the end.",
        "We're finished. All set. Bye. Done",
    ]
    result = await aggregate_chunked_results(func, chunked, kernel, KernelArguments())
    print(result)
    assert result == "\n".join(chunked)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
