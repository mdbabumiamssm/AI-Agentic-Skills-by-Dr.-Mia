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

# Prerequisites

Before you begin, ensure you have the following prerequisites:

- .NET 10 SDK or later
- Docker installed and running on your machine
- An Ollama model downloaded into Ollama

To download and start Ollama on Docker using CPU, run the following command in your terminal.

```powershell
docker run -d -v "c:\temp\ollama:/root/.ollama" -p 11434:11434 --name ollama ollama/ollama
```

To download and start Ollama on Docker using GPU, run the following command in your terminal.

```powershell
docker run -d --gpus=all -v "c:\temp\ollama:/root/.ollama" -p 11434:11434 --name ollama ollama/ollama
```

After the container has started, launch a Terminal window for the docker container, e.g. if using docker desktop, choose Open in Terminal from actions.

From this terminal download the required models, e.g. here we are downloading the phi3 model.

```text
ollama pull gpt-oss
```

Set the following environment variables:

```powershell
$env:OLLAMA_ENDPOINT="http://localhost:11434"
$env:OLLAMA_MODEL_NAME="gpt-oss"
```


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->