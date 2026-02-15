# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""External preprocessing functions."""

from __future__ import annotations

from sklearn.utils import deprecated

from ...preprocessing import _scrublet
from ._bbknn import bbknn
from ._dca import dca
from ._harmony_integrate import harmony_integrate
from ._hashsolo import hashsolo
from ._magic import magic
from ._mnn_correct import mnn_correct
from ._scanorama_integrate import scanorama_integrate

scrublet = deprecated("Import from sc.pp instead")(_scrublet.scrublet)
scrublet_simulate_doublets = deprecated("Import from sc.pp instead")(
    _scrublet.scrublet_simulate_doublets
)

__all__ = [
    "bbknn",
    "dca",
    "harmony_integrate",
    "hashsolo",
    "magic",
    "mnn_correct",
    "scanorama_integrate",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
