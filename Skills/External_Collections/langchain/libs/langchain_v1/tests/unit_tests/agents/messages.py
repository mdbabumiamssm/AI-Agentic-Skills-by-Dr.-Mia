# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Redefined messages as a work-around for pydantic issue with AnyStr.

The code below creates version of pydantic models
that will work in unit tests with AnyStr as id field
Please note that the `id` field is assigned AFTER the model is created
to workaround an issue with pydantic ignoring the __eq__ method on
subclassed strings.
"""

from typing import Any

from langchain_core.messages import HumanMessage, ToolMessage

from .any_str import AnyStr


def _AnyIdHumanMessage(**kwargs: Any) -> HumanMessage:  # noqa: N802
    """Create a human message with an any id field."""
    message = HumanMessage(**kwargs)
    message.id = AnyStr()
    return message


def _AnyIdToolMessage(**kwargs: Any) -> ToolMessage:  # noqa: N802
    """Create a tool message with an any id field."""
    message = ToolMessage(**kwargs)
    message.id = AnyStr()
    return message

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
