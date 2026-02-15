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

# Starting With Semantic Kernel Vector Stores

This project contains a step by step guide to get started using Vector Stores with the Semantic Kernel.

The examples can be run as integration tests but their code can also be copied to stand-alone programs.

## Configuring Secrets

Most of the examples will require secrets and credentials, to access OpenAI, Azure OpenAI,
Vector Stores and other resources. We suggest using .NET
[Secret Manager](https://learn.microsoft.com/aspnet/core/security/app-secrets)
to avoid the risk of leaking secrets into the repository, branches and pull requests.
You can also use environment variables if you prefer.

To set your secrets with Secret Manager:

```
cd dotnet/samples/GettingStartedWithVectorStores

dotnet user-secrets init

dotnet user-secrets set "AzureOpenAIEmbeddings:DeploymentName" "..."
dotnet user-secrets set "AzureOpenAIEmbeddings:Endpoint" "..."

dotnet user-secrets set "AzureAISearch:Endpoint" "..."
dotnet user-secrets set "AzureAISearch:ApiKey" "..."
```

To set your secrets with environment variables, use these names:

```
AzureOpenAIEmbeddings__DeploymentName
AzureOpenAIEmbeddings__Endpoint

AzureAISearch__Endpoint
AzureAISearch__ApiKey
```


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->