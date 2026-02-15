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
from typing import Any, ClassVar

import boto3
from botocore.config import Config

from semantic_kernel.connectors.ai.bedrock.services.model_provider.bedrock_model_provider import BedrockModelProvider
from semantic_kernel.kernel_pydantic import KernelBaseModel


class BedrockBase(KernelBaseModel, ABC):
    """Amazon Bedrock Service Base Class."""

    MODEL_PROVIDER_NAME: ClassVar[str] = "bedrock"

    # Amazon Bedrock Clients
    # Runtime Client: Use for inference
    bedrock_runtime_client: Any
    # Client: Use for model management
    bedrock_client: Any

    bedrock_model_provider: BedrockModelProvider | None = None

    def __init__(
        self,
        *,
        runtime_client: Any | None = None,
        client: Any | None = None,
        bedrock_model_provider: BedrockModelProvider | None = None,
        **kwargs: Any,
    ) -> None:
        """Initialize the Amazon Bedrock Base Class.

        Args:
            runtime_client: The Amazon Bedrock runtime client to use.
            client: The Amazon Bedrock client to use.
            bedrock_model_provider: The Bedrock model provider to use.
                If not provided, the model provider will be extracted from the model ID.
                When using an Application Inference Profile where the model provider is not part
                of the model ID, this setting must be provided.
            **kwargs: Additional keyword arguments.
        """
        config = Config(user_agent_extra="x-client-framework:semantic-kernel")

        super().__init__(
            bedrock_runtime_client=runtime_client or boto3.client("bedrock-runtime", config=config),
            bedrock_client=client or boto3.client("bedrock"),
            bedrock_model_provider=bedrock_model_provider,
            **kwargs,
        )

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
