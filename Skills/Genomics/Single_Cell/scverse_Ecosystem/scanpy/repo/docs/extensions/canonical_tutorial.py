# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Extension for a stub ``canonical-tutorial`` directive."""

from __future__ import annotations

from typing import TYPE_CHECKING

from sphinx.util.docutils import SphinxDirective

if TYPE_CHECKING:
    from typing import ClassVar

    from docutils import nodes
    from sphinx.application import Sphinx


class CanonicalTutorial(SphinxDirective):
    """In the scanpy-tutorials repo, this links to the canonical location (here!)."""

    required_arguments: ClassVar = 1

    def run(self) -> list[nodes.Node]:  # noqa: D102
        return []


def setup(app: Sphinx) -> None:
    """App setup hook."""
    app.add_directive("canonical-tutorial", CanonicalTutorial)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
