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

from dotenv import dotenv_values
from promptflow import PFClient
from promptflow.entities import AzureOpenAIConnection

pf_client = PFClient()

# Run a single test of a flow
#################################################

# Load the configuration from the .env file
config = dotenv_values(".env")
deployment_type = config.get("AZURE_OPENAI_DEPLOYMENT_TYPE", None)
if deployment_type == "chat-completion":
    deployment_name = config.get("AZURE_OPENAI_CHAT_COMPLETION_DEPLOYMENT_NAME", None)
elif deployment_type == "text-completion":
    deployment_name = config.get("AZURE_OPENAI_TEXT_COMPLETION_DEPLOYMENT_NAME", None)

# Define the inputs of the flow
inputs = {
    "text": "What is 2 plus 3?",
    "deployment_type": deployment_type,
    "deployment_name": deployment_name,
}

# Initialize an AzureOpenAIConnection object
connection = AzureOpenAIConnection(
    name="AzureOpenAIConnection",
    type="Custom",
    api_key=config.get("AZURE_OPENAI_API_KEY", None),
    api_base=config.get("AZURE_OPENAI_ENDPOINT", None),
    api_version="2023-03-15-preview",
)

# Add connections to the Prompt flow client
pf_client.connections.create_or_update(connection)

# Run the flow
flow_result = pf_client.test(flow="perform_math", inputs=inputs)

# Print the outputs of the flow
print(f"Flow outputs: {flow_result}")

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
