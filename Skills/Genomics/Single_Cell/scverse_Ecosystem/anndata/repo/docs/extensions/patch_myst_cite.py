# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Override MyST’s cite role with one that works."""

from __future__ import annotations

from types import MappingProxyType
from typing import TYPE_CHECKING

from docutils import nodes, utils

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence
    from typing import Any

    from docutils.parsers.rst.states import Inliner
    from sphinx.application import Sphinx


def cite_role(  # noqa: PLR0917
    name: str,
    rawsource: str,
    text: str,
    lineno: int,
    inliner: Inliner,
    options: Mapping[str, Any] = MappingProxyType({}),
    content: Sequence[str] = (),
) -> tuple[list[nodes.Node], list[nodes.system_message]]:
    key = utils.unescape(text)
    node = nodes.citation_reference(f"[{key}]_", key)
    return [node], []


def setup(app: Sphinx):
    app.add_role("cite", cite_role, override=True)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
