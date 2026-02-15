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

doc_neighbors_key = """\
neighbors_key
    Where to look for neighbors connectivities.
    If not specified, this retrieves ``.obsp['connectivities']`` for connectivities
    (default storage place for :func:`~scanpy.pp.neighbors`).
    If specified, this retrieves
    ``.obsp[.uns[neighbors_key]['connectivities_key']]`` for connectivities.
"""

doc_use_rep = """\
use_rep
    Use the indicated representation. `'X'` or any key for `.obsm` is valid.
    If `None`, the representation is chosen automatically:
    For `.n_vars` < :attr:`~scanpy.settings.N_PCS` (default: 50), `.X` is used, otherwise 'X_pca' is used.
    If 'X_pca' is not present, it’s computed with default parameters or `n_pcs` if present.\
"""

doc_n_pcs = """\
n_pcs
    Use this many PCs. If `n_pcs==0` use `.X` if `use_rep is None`.\
"""

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
