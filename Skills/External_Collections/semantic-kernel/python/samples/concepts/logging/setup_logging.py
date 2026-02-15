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

import logging
import os

from semantic_kernel import Kernel
from semantic_kernel.connectors.ai.open_ai import OpenAIChatCompletion
from semantic_kernel.utils.logging import setup_logging


async def main():
    setup_logging()

    # Set the logging level for  semantic_kernel.kernel to DEBUG.
    logging.getLogger("kernel").setLevel(logging.DEBUG)

    kernel = Kernel()

    service_id = "chat-gpt"
    kernel.add_service(OpenAIChatCompletion(service_id=service_id, ai_model_id="gpt-3.5-turbo"))

    plugins_directory = os.path.join(__file__, "../../../../../prompt_template_samples/")
    plugin = kernel.add_plugin(parent_directory=plugins_directory, plugin_name="FunPlugin")

    joke_function = plugin["Joke"]

    result = await kernel.invoke(joke_function, input="time travel to dinosaur age", style="silly")

    print(result)


if __name__ == "__main__":
    import asyncio

    asyncio.run(main())

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
