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

﻿
# A2A Client Sample
Show how to create an A2A Client with a command line interface which invokes agents using the A2A protocol.

## Run the Sample

To run the sample, follow these steps:

1. Run the A2A client:
    ```bash
    cd A2AClient
    dotnet run
    ```  
2. Enter your request e.g. "Show me all invoices for Contoso?"

## Set Secrets with Secret Manager

The agent urls are provided as a ` ` delimited list of strings

```text
cd dotnet/samples/Demos/A2AClientServer/A2AClient

dotnet user-secrets set "A2AClient:ModelId" "..."
dotnet user-secrets set "A2AClient":ApiKey" "..."
dotnet user-secrets set "A2AClient:AgentUrls" "http://localhost:5000/policy;http://localhost:5000/invoice;http://localhost:5000/logistics"
```

<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->