<!--
# COPYRIGHT NOTICE
# This file is part of the "Universal AI Agentic Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

-->

# Sample Integration Tests

Integration tests that validate the Durable Agent Framework samples by running them as Azure Functions.

## Setup

### 1. Create `.env` file

Copy `.env.example` to `.env` and fill in your Azure credentials:

```bash
cp .env.example .env
```

Required variables:
- `AZURE_OPENAI_ENDPOINT`
- `AZURE_OPENAI_CHAT_DEPLOYMENT_NAME`
- `AZURE_OPENAI_API_KEY`
- `AzureWebJobsStorage`
- `DURABLE_TASK_SCHEDULER_CONNECTION_STRING`
- `FUNCTIONS_WORKER_RUNTIME`

### 2. Start required services

**Azurite (for orchestration tests):**
```bash
docker run -d -p 10000:10000 -p 10001:10001 -p 10002:10002 mcr.microsoft.com/azure-storage/azurite
```

**Durable Task Scheduler:**
```bash
docker run -d -p 8080:8080 -p 8082:8082 -e DTS_USE_DYNAMIC_TASK_HUBS=true mcr.microsoft.com/dts/dts-emulator:latest
```

## Running Tests

The tests automatically start and stop the Azure Functions app for each sample.

### Run all sample tests
```bash
uv run pytest packages/azurefunctions/tests/integration_tests -v
```

### Run specific sample
```bash
uv run pytest packages/azurefunctions/tests/integration_tests/test_01_single_agent.py -v
```

### Run with verbose output
```bash
uv run pytest packages/azurefunctions/tests/integration_tests -sv
```

## How It Works

Each test file uses pytest markers to automatically configure and start the function app:

```python
pytestmark = [
    pytest.mark.sample("01_single_agent"),
    pytest.mark.usefixtures("function_app_for_test"),
    skip_if_azure_functions_integration_tests_disabled,
]
```

The `function_app_for_test` fixture:
1. Loads environment variables from `.env`
2. Validates required variables are present
3. Starts the function app on a dynamically allocated port
4. Waits for the app to be ready
5. Runs your tests
6. Tears down the function app

## Troubleshooting


**Missing environment variables:**
Ensure your `.env` file contains all required variables from `.env.example`.

**Tests timeout:**
Check that Azure OpenAI credentials are valid and the service is accessible.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->