# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Extension to skip inherited methods and properties in autosummary."""

from __future__ import annotations

from traceback import walk_stack
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Literal

    from sphinx.application import Sphinx
    from sphinx.ext.autodoc import Options


def skip_inherited(  # noqa: PLR0917
    app: Sphinx,
    what: Literal[
        "module", "class", "exception", "function", "method", "attribute", "property"
    ],
    name: str,
    obj: object,
    skip: bool,  # noqa: FBT001
    options: Options | dict[str, object],
) -> bool | None:
    """Skip inherited members."""
    # Skip `getdoc` property
    if what == "method" and name == "getdoc":
        return True

    # find parent class
    for frame, _ in walk_stack(None):
        # Sadly `autodoc-skip-member` doesn’t give access to the parent object, so we need to do this.
        # Find this stack frame: https://github.com/sphinx-doc/sphinx/blob/a5366394ae527712c4edfeb07a5fbeecd4ca72e1/sphinx/ext/autosummary/generate.py#L496-L517
        if frame.f_code.co_name == "_get_members" and frame.f_code.co_filename.endswith(
            "/generate.py"
        ):
            parent = frame.f_locals["obj"]
            if not isinstance(parent, type):
                return None
            break
    else:
        return None

    # return if it’s a member of the parent class
    typ = parent
    while typ is not type:
        if name in typ.__dict__:
            return None
        # Metaclasses need the while loop to trace back to their parents.
        # See `SingletonMeta` and the documentation of `settings` for an example.
        typ = type(typ)

    # Skip this `name` because the `parent` `type` lacks the `name`, which indicates `name` does not belong to `parent`
    return True


def setup(app: Sphinx) -> None:
    """App setup hook."""
    app.connect("autodoc-skip-member", skip_inherited)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
