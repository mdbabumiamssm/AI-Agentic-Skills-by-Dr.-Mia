# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from .astronauts import get_num_astronauts
from .astronauts import get_coords_iss

def test_astro():
    assert type(get_num_astronauts())==int

def test_iss():
    latitude, longitude = get_coords_iss()
    assert type(latitude)==float
    assert type(longitude)==float

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
