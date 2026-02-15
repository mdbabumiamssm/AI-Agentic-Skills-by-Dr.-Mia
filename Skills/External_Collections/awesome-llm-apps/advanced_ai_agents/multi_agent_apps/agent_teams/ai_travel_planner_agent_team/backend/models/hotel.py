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
from typing import List

class HotelResult(BaseModel):
    hotel_name: str = Field(description="The name of the hotel")
    price: str = Field(description="The price of the hotel")
    rating: str = Field(description="The rating of the hotel")
    address: str = Field(description="The address of the hotel")
    amenities: List[str] = Field(description="The amenities of the hotel")
    description: str = Field(description="The description of the hotel")
    url: str = Field(description="The url of the hotel")

class HotelResults(BaseModel):
    hotels: List[HotelResult] = Field(description="The list of hotels")

class HotelSearchRequest(BaseModel):
    destination: str = Field(description="The destination city or area")
    check_in: str = Field(description="The date of check-in in the format 'YYYY-MM-DD'")
    check_out: str = Field(description="The date of check-out in the format 'YYYY-MM-DD'")
    adults: int = Field(description="The number of adults")
    children: int = Field(description="The number of children")
    rooms: int = Field(description="The number of rooms")
    sort: str = Field(description="The sort order")
__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
