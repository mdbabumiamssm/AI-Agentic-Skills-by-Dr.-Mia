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


from semantic_kernel.functions.kernel_arguments import KernelArguments
from semantic_kernel.functions.kernel_function import KernelFunction
from semantic_kernel.kernel import Kernel


async def aggregate_chunked_results(
    func: KernelFunction, chunked_results: list[str], kernel: Kernel, arguments: KernelArguments
) -> str:
    """Aggregate the results from the chunked results."""
    results = []
    for chunk in chunked_results:
        arguments["input"] = chunk
        result = await func.invoke(kernel, arguments)

        results.append(str(result))

    return "\n".join(results)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
