# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Models for the Constitutional AI chain."""

from pydantic import BaseModel


class ConstitutionalPrinciple(BaseModel):
    """Class for a constitutional principle."""

    critique_request: str
    revision_request: str
    name: str = "Constitutional Principle"

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
