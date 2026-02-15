# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from pydantic import BaseModel, Field


class Pagination(BaseModel):
    total_items: int = Field(description="Total number of items.", examples=[42])
    total_pages: int = Field(description="Total number of pages.", examples=[97])
    current_page: int = Field(description="Current_page page number.", examples=[1])
    page_size: int = Field(description="Number of items per page.", examples=[25])

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
