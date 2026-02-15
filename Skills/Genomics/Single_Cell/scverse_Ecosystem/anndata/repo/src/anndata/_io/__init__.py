# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from __future__ import annotations

from ..utils import warn

__all__: list[str] = []


def __getattr__(key: str):
    from .. import io

    attr = getattr(io, key)
    msg = (
        f"Importing {key} from `anndata._io` is deprecated. "
        "Please use anndata.io instead."
    )
    warn(msg, FutureWarning)
    return attr

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
