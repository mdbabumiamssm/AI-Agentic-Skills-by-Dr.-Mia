# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from agno.agent import Agent
from agno.models.ollama import Ollama
from agno.playground import Playground, serve_playground_app

reasoning_agent = Agent(name="Reasoning Agent", model=Ollama(id="qwq:32b"), markdown=True)

# UI for Reasoning agent
app = Playground(agents=[reasoning_agent]).get_app()

# Run the Playground app
if __name__ == "__main__":
    serve_playground_app("local_ai_reasoning_agent:app", reload=True)
__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
