# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from langchain_core.utils.pydantic import PYDANTIC_VERSION


def get_pydantic_major_version() -> int:
    """Get the major version of Pydantic.

    Returns:
        The major version of Pydantic.
    """
    return PYDANTIC_VERSION.major


__all__ = ["get_pydantic_major_version"]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
