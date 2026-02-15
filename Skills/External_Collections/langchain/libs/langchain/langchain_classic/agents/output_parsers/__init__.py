# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Parsing utils to go from string to AgentAction or Agent Finish.

AgentAction means that an action should be taken.
This contains the name of the tool to use, the input to pass to that tool,
and a `log` variable (which contains a log of the agent's thinking).

AgentFinish means that a response should be given.
This contains a `return_values` dictionary. This usually contains a
single `output` key, but can be extended to contain more.
This also contains a `log` variable (which contains a log of the agent's thinking).
"""

from langchain_classic.agents.output_parsers.json import JSONAgentOutputParser
from langchain_classic.agents.output_parsers.openai_functions import (
    OpenAIFunctionsAgentOutputParser,
)
from langchain_classic.agents.output_parsers.react_json_single_input import (
    ReActJsonSingleInputOutputParser,
)
from langchain_classic.agents.output_parsers.react_single_input import (
    ReActSingleInputOutputParser,
)
from langchain_classic.agents.output_parsers.self_ask import SelfAskOutputParser
from langchain_classic.agents.output_parsers.tools import ToolsAgentOutputParser
from langchain_classic.agents.output_parsers.xml import XMLAgentOutputParser

__all__ = [
    "JSONAgentOutputParser",
    "OpenAIFunctionsAgentOutputParser",
    "ReActJsonSingleInputOutputParser",
    "ReActSingleInputOutputParser",
    "SelfAskOutputParser",
    "ToolsAgentOutputParser",
    "XMLAgentOutputParser",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
