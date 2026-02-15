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

from collections.abc import Awaitable, Callable
from typing import Any
from urllib.parse import urlparse

import httpx
from pydantic import Field

from semantic_kernel.connectors.openapi_plugin.operation_selection_predicate_context import (
    OperationSelectionPredicateContext,
)
from semantic_kernel.kernel_pydantic import KernelBaseModel

AuthCallbackType = Callable[..., Awaitable[Any]]


class OpenAPIFunctionExecutionParameters(KernelBaseModel):
    """OpenAPI function execution parameters."""

    http_client: httpx.AsyncClient | None = None
    auth_callback: AuthCallbackType | None = None
    server_url_override: str | None = None
    ignore_non_compliant_errors: bool = False
    user_agent: str | None = None
    enable_dynamic_payload: bool = True
    enable_payload_namespacing: bool = False
    operations_to_exclude: list[str] = Field(default_factory=list, description="The operationId(s) to exclude")
    operation_selection_predicate: Callable[[OperationSelectionPredicateContext], bool] | None = None
    timeout: float | None = Field(
        None, description="Default timeout in seconds for HTTP requests. Uses httpx default (5 seconds) if None."
    )

    def model_post_init(self, __context: Any) -> None:
        """Post initialization method for the model."""
        from semantic_kernel.utils.telemetry.user_agent import HTTP_USER_AGENT

        if self.server_url_override:
            parsed_url = urlparse(self.server_url_override)
            if not parsed_url.scheme or not parsed_url.netloc:
                raise ValueError(f"Invalid server_url_override: {self.server_url_override}")

        if not self.user_agent:
            self.user_agent = HTTP_USER_AGENT

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
