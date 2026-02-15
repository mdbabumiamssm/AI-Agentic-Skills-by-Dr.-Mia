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

from semantic_kernel.connectors.ai.prompt_execution_settings import PromptExecutionSettings
from semantic_kernel.functions.function_result import FunctionResult
from semantic_kernel.kernel_pydantic import KernelBaseModel
from semantic_kernel.services.ai_service_client_base import AIServiceClientBase


class PromptRenderingResult(KernelBaseModel):
    """Represents the result of rendering a prompt template.

    Attributes:
        rendered_prompt (str): The rendered prompt.
        ai_service (Any): The AI service that rendered the prompt.
        execution_settings (PromptExecutionSettings): The execution settings for the prompt.
        function_result (FunctionResult): The result of executing the prompt.
    """

    rendered_prompt: str
    ai_service: AIServiceClientBase
    execution_settings: PromptExecutionSettings
    function_result: FunctionResult | None = None

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
