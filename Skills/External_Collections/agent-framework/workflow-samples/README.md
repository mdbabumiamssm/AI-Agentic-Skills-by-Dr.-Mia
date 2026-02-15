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

# Declarative Workflows

A _Declarative Workflow_ is defined as a single YAML file and
may be executed locally no different from any regular `Workflow` that is defined by code.

The difference is that the workflow definition is loaded from a YAML file instead of being defined in code:

```c#
Workflow workflow = DeclarativeWorkflowBuilder.Build("Marketing.yaml", options);
```

These example workflows may be executed by the workflow
[Samples](../dotnet/samples/GettingStarted/Workflows/Declarative)
that are present in this repository.

> See the [README.md](../dotnet/samples/GettingStarted/Workflows/Declarative/README.md) 
 associated with the samples for configuration details.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->