# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Different methods for rendering Tools to be passed to LLMs.

Depending on the LLM you are using and the prompting strategy you are using,
you may want Tools to be rendered in a different way.
This module contains various ways to render tools.
"""

# For backwards compatibility
from langchain_core.tools import (
    render_text_description,
    render_text_description_and_args,
)
from langchain_core.utils.function_calling import (
    convert_to_openai_function as format_tool_to_openai_function,
)
from langchain_core.utils.function_calling import (
    convert_to_openai_tool as format_tool_to_openai_tool,
)

__all__ = [
    "format_tool_to_openai_function",
    "format_tool_to_openai_tool",
    "render_text_description",
    "render_text_description_and_args",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
