# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from langchain_core.agents import AgentAction

from langchain_classic.agents.conversational.output_parser import ConvoOutputParser


def test_normal_output_parsing() -> None:
    _test_convo_output(
        """
Action: my_action
Action Input: my action input
""",
        "my_action",
        "my action input",
    )


def test_multiline_output_parsing() -> None:
    _test_convo_output(
        """
Thought: Do I need to use a tool? Yes
Action: evaluate_code
Action Input: Evaluate Code with the following Python content:
```python
print("Hello fifty shades of gray mans!"[::-1])  # noqa: T201
```
""",
        "evaluate_code",
        """
Evaluate Code with the following Python content:
```python
print("Hello fifty shades of gray mans!"[::-1])  # noqa: T201
```""".lstrip(),
    )


def _test_convo_output(text: str, expected_tool: str, expected_tool_input: str) -> None:
    result = ConvoOutputParser().parse(text.strip())
    assert isinstance(result, AgentAction)
    assert result.tool == expected_tool
    assert result.tool_input == expected_tool_input

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
