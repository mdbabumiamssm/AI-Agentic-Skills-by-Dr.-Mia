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
import os

from semantic_kernel import Kernel
from semantic_kernel.functions.kernel_arguments import KernelArguments


async def main():
    """OpenAPI Sample Client"""
    kernel = Kernel()

    spec_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))),
        "plugins",
        "openapi",
        "openapi.yaml",
    )

    openapi_plugin = kernel.add_plugin_from_openapi(plugin_name="openApiPlugin", openapi_document_path=spec_path)

    arguments = KernelArguments(
        input="hello world",
        name="John",
        q=0.7,
        Header="example-header",
    )

    result = await kernel.invoke(openapi_plugin["helloWorld"], arguments=arguments)

    print(result)


if __name__ == "__main__":
    asyncio.run(main())

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
