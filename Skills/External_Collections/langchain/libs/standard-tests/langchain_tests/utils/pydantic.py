# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Utilities for working with pydantic models."""


def get_pydantic_major_version() -> int:
    """Get the major version of Pydantic."""
    try:
        import pydantic  # noqa: PLC0415

        return int(pydantic.__version__.split(".")[0])
    except ImportError:
        return 0


PYDANTIC_MAJOR_VERSION = get_pydantic_major_version()

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
