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

if TYPE_CHECKING:
    from torch import Tensor


def rowwise_max_excluding_diagonal(matrix: Tensor) -> Tensor:
    """Get the rowwise maximum of a matrix excluding the diagonal."""
    import torch

    assert matrix.ndim == 2
    num_cols = matrix.shape[1]
    mask = (1 - torch.eye(num_cols, device=matrix.device)).bool()
    return (
        torch.where(mask, matrix, torch.tensor(-float("inf"), device=matrix.device))
        .max(axis=1)
        .values
    )


def simple_reciprocal(w: Tensor, eps: float = 1e-6) -> Tensor:
    """Convert distances to similarities via a reciprocal."""
    return 1.0 / (w + eps)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
