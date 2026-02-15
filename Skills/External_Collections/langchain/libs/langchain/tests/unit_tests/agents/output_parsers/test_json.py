# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from langchain_core.agents import AgentAction, AgentFinish

from langchain_classic.agents.output_parsers.json import JSONAgentOutputParser


def test_tool_usage() -> None:
    parser = JSONAgentOutputParser()
    _input = """    ```
{
  "action": "search",
  "action_input": "2+2"
}
```"""
    output = parser.invoke(_input)
    expected_output = AgentAction(tool="search", tool_input="2+2", log=_input)
    assert output == expected_output


def test_finish() -> None:
    parser = JSONAgentOutputParser()
    _input = """```
{
  "action": "Final Answer",
  "action_input": "4"
}
```"""
    output = parser.invoke(_input)
    expected_output = AgentFinish(return_values={"output": "4"}, log=_input)
    assert output == expected_output

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
