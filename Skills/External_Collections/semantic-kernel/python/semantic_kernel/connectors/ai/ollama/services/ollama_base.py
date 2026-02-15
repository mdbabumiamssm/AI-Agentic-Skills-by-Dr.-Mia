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

from abc import ABC
from typing import ClassVar

from ollama import AsyncClient

from semantic_kernel.kernel_pydantic import KernelBaseModel


class OllamaBase(KernelBaseModel, ABC):
    """Ollama service base.

    Args:
        client [AsyncClient]: An Ollama client to use for the service.
    """

    MODEL_PROVIDER_NAME: ClassVar[str] = "ollama"

    client: AsyncClient

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
