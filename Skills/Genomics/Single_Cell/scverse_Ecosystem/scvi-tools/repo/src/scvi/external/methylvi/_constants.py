# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from typing import NamedTuple


class _METHYLVI_REGISTRY_KEYS_NT(NamedTuple):
    MC_KEY: str = "mc"
    COV_KEY: str = "cov"


METHYLVI_REGISTRY_KEYS = _METHYLVI_REGISTRY_KEYS_NT()

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
