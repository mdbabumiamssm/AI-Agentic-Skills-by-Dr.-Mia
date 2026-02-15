# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Extension to skip deprecated methods and properties in autosummary."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Literal

    from sphinx.application import Sphinx
    from sphinx.ext.autodoc import Options


def skip_deprecated(  # noqa: PLR0917
    app: Sphinx,
    what: Literal[
        "module", "class", "exception", "function", "method", "attribute", "property"
    ],
    name: str,
    obj: object,
    skip: bool,  # noqa: FBT001
    options: Options | dict[str, object],
) -> bool | None:
    """Skip deprecated members."""
    if hasattr(obj, "__deprecated__"):
        return True
    return None


def setup(app: Sphinx) -> None:
    """App setup hook."""
    app.connect("autodoc-skip-member", skip_deprecated)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
