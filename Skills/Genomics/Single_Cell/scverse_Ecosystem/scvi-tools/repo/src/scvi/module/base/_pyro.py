# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from pyro.infer import Predictive

from ._decorators import auto_move_data


class AutoMoveDataPredictive(Predictive):
    """Auto move data."""

    @auto_move_data
    def forward(self, *args, **kwargs):
        return super().forward(*args, **kwargs)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
