# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from typing import Iterator
from agno.agent import Agent, RunResponse
from agno.models.openai import OpenAIChat
from agno.utils.pprint import pprint_run_response
from dotenv import load_dotenv

load_dotenv()
agent = Agent(model=OpenAIChat(id="gpt-4o-mini"))
response: RunResponse = agent.run("Tell me a 5 second short story about a robot")
response_stream: Iterator[RunResponse] = agent.run("Tell me a 5 second short story about a lion", stream=True)
pprint_run_response(response, markdown=True)
pprint_run_response(response_stream, markdown=True)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
