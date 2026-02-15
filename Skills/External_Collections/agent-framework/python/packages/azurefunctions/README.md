<!--
# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

-->

# Get Started with Microsoft Agent Framework Durable Functions

[![PyPI](https://img.shields.io/pypi/v/agent-framework-azurefunctions)](https://pypi.org/project/agent-framework-azurefunctions/)

Please install this package via pip:

```bash
pip install agent-framework-azurefunctions --pre
```

## Durable Agent Extension

The durable agent extension lets you host Microsoft Agent Framework agents on Azure Durable Functions so they can persist state, replay conversation history, and recover from failures automatically.

### Basic Usage Example

See the durable functions integration sample in the repository to learn how to:

```python
from agent_framework.azure import AgentFunctionApp

_app = AgentFunctionApp()
```

- Register agents with `AgentFunctionApp`
- Post messages using the generated `/api/agents/{agent_name}/run` endpoint

For more details, review the Python [README](https://github.com/microsoft/agent-framework/tree/main/python/README.md) and the samples directory.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->