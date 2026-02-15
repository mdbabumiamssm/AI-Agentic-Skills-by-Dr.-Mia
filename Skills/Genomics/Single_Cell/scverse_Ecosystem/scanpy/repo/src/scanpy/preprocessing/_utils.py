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

from typing import TYPE_CHECKING

import numpy as np
from sklearn.random_projection import sample_without_replacement

if TYPE_CHECKING:
    from typing import Literal

    from numpy.typing import NDArray

    from .._utils.random import _LegacyRandom


def sample_comb(
    dims: tuple[int, ...],
    nsamp: int,
    *,
    random_state: _LegacyRandom = None,
    method: Literal[
        "auto", "tracking_selection", "reservoir_sampling", "pool"
    ] = "auto",
) -> NDArray[np.int64]:
    """Randomly sample indices from a grid, without repeating the same tuple."""
    idx = sample_without_replacement(
        np.prod(dims), nsamp, random_state=random_state, method=method
    )
    return np.vstack(np.unravel_index(idx, dims)).T

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
