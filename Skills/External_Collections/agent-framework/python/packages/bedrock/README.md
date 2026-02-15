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

# Get Started with Microsoft Agent Framework Bedrock

Install the provider package:

```bash
pip install agent-framework-bedrock --pre
```

## Bedrock Integration

The Bedrock integration enables Microsoft Agent Framework applications to call Amazon Bedrock models with familiar chat abstractions, including tool/function calling when you attach tools through `ChatOptions`.

### Basic Usage Example

See the [Bedrock sample script](samples/bedrock_sample.py) for a runnable end-to-end script that:

- Loads credentials from the `BEDROCK_*` environment variables
- Instantiates `BedrockChatClient`
- Sends a simple conversation turn and prints the response


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->