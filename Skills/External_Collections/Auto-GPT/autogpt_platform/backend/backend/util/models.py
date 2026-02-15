# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""
Shared models and types used across the backend to avoid circular imports.
"""

import pydantic


class Pagination(pydantic.BaseModel):
    total_items: int = pydantic.Field(
        description="Total number of items.", examples=[42]
    )
    total_pages: int = pydantic.Field(
        description="Total number of pages.", examples=[2]
    )
    current_page: int = pydantic.Field(
        description="Current_page page number.", examples=[1]
    )
    page_size: int = pydantic.Field(
        description="Number of items per page.", examples=[25]
    )

    @staticmethod
    def empty() -> "Pagination":
        return Pagination(
            total_items=0,
            total_pages=0,
            current_page=0,
            page_size=0,
        )

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
