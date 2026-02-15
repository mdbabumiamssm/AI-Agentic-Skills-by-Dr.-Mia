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

"""This module defines an enumeration representing different services."""

from enum import Enum


class Service(Enum):
    """Attributes:
    OpenAI (str): Represents the OpenAI service.
    AzureOpenAI (str): Represents the Azure OpenAI service.
    HuggingFace (str): Represents the HuggingFace service.
    """

    OpenAI = "openai"
    AzureOpenAI = "azureopenai"
    HuggingFace = "huggingface"

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
