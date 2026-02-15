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

from typing_extensions import deprecated

from semantic_kernel.connectors.ai.google.vertex_ai.vertex_ai_settings import VertexAISettings
from semantic_kernel.kernel_pydantic import KernelBaseModel


@deprecated("VertexAIBase is deprecated and will be removed after 01/01/2026. Use google_ai connectors instead.")
class VertexAIBase(KernelBaseModel, ABC):
    """Vertex AI Service."""

    MODEL_PROVIDER_NAME: ClassVar[str] = "vertexai"

    service_settings: VertexAISettings

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
