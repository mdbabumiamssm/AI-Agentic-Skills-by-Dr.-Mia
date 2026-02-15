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

### Running the OpenAPI syntax example

For more generic setup instructions, including how to install the `uv` tool, see the [main README](../../../../DEV_SETUP.md).

1. In a terminal, navigate to `semantic_kernel/python/samples/concepts/plugins/openapi`.

2. Run `uv sync` followed by `source .venv/bin/activate` to enter the virtual environment (depending on the os, the activate script may be in a different location).

3. Start the server by running `python openapi_server.py`.

4. In another terminal, do steps 1 & 2. Then, run `python openapi_client.py`, which will register a plugin representing the API defined in openapi.yaml


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->