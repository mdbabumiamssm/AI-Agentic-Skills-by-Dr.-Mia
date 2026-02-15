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

"""Backend tool rendering endpoint."""

from agent_framework.ag_ui import add_agent_framework_fastapi_endpoint
from agent_framework.azure import AzureOpenAIChatClient
from fastapi import FastAPI

from ...agents.weather_agent import weather_agent


def register_backend_tool_rendering(app: FastAPI) -> None:
    """Register the backend tool rendering endpoint.

    Args:
        app: The FastAPI application.
    """
    # Create a chat client and call the factory function
    chat_client = AzureOpenAIChatClient()

    add_agent_framework_fastapi_endpoint(
        app,
        weather_agent(chat_client),
        "/backend_tool_rendering",
    )

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
