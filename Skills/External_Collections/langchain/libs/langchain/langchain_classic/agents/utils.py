# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from collections.abc import Sequence

from langchain_core.tools import BaseTool


def validate_tools_single_input(class_name: str, tools: Sequence[BaseTool]) -> None:
    """Validate tools for single input.

    Args:
        class_name: Name of the class.
        tools: List of tools to validate.

    Raises:
        ValueError: If a multi-input tool is found in tools.
    """
    for tool in tools:
        if not tool.is_single_input:
            msg = f"{class_name} does not support multi-input tool {tool.name}."
            raise ValueError(msg)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
