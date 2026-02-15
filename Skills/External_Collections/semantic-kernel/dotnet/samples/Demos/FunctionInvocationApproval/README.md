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

# Function Invocation Approval

This console application shows how to use function invocation filter (`IFunctionInvocationFilter`) to invoke a Kernel Function only if such operation was approved.
If function invocation was rejected, the result will contain the reason why, so the LLM can respond appropriately.

The application uses a sample plugin which builds software by following these development stages: collection of requirements, design, implementation, testing and deployment.

Each step can be approved or rejected. Based on that, the LLM will decide how to proceed.

## Configuring Secrets

The example requires credentials to access OpenAI or Azure OpenAI.

If you have set up those credentials as secrets within Secret Manager or through environment variables for other samples from the solution in which this project is found, they will be re-used.

### To set your secrets with Secret Manager:

```
cd dotnet/samples/Demos/FunctionInvocationApproval

dotnet user-secrets init

dotnet user-secrets set "OpenAI:ChatModelId" "..."
dotnet user-secrets set "OpenAI:ApiKey" "..."

dotnet user-secrets set "AzureOpenAI:ChatDeploymentName" "..."
dotnet user-secrets set "AzureOpenAI:Endpoint" "https://... .openai.azure.com/"
dotnet user-secrets set "AzureOpenAI:ApiKey" "..."
```

### To set your secrets with environment variables

Use these names:

```
# OpenAI
OpenAI__ChatModelId
OpenAI__ApiKey

# Azure OpenAI
AzureOpenAI__ChatDeploymentName
AzureOpenAI__Endpoint
AzureOpenAI__ApiKey
```


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->