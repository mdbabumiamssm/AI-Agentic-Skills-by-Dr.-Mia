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

from typing import Literal

from semantic_kernel.kernel_pydantic import KernelBaseSettings


class ServiceSettings(KernelBaseSettings):
    """The Learn Resources Service Settings.

    The settings are first loaded from environment variables. If the
    environment variables are not found, the settings can be loaded from a .env file with the
    encoding 'utf-8' as default or the specific encoding. If the settings are not found in the
    .env file, the settings are ignored; however, validation will fail alerting that the settings
    are missing.

    Args:
        global_llm_service: The LLM service to use for the samples, either "OpenAI" or "AzureOpenAI"
            If not provided, defaults to "AzureOpenAI".
    """

    global_llm_service: Literal["OpenAI", "AzureOpenAI"] = "AzureOpenAI"

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
