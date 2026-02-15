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

# Prerequisites

WARNING: ONNX doesn't support function calling, so any function tools passed to the agent will be ignored.

Before you begin, ensure you have the following prerequisites:

- .NET 10 SDK or later
- An ONNX model downloaded to your machine

You can download an ONNX model from hugging face, using git clone:

```powershell
git clone https://huggingface.co/microsoft/Phi-4-mini-instruct-onnx
```

Set the following environment variables:

```powershell
$env:ONNX_MODEL_PATH="C:\repos\Phi-4-mini-instruct-onnx\cpu_and_mobile\cpu-int4-rtn-block-32-acc-level-4" # Replace with your model path
```


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->