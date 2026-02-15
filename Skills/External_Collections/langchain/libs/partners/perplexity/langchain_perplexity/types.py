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

from typing import Literal

from pydantic import BaseModel


class UserLocation(BaseModel):
    latitude: float | None = None
    longitude: float | None = None
    country: str | None = None
    region: str | None = None
    city: str | None = None


class WebSearchOptions(BaseModel):
    search_context_size: Literal["low", "medium", "high"] | None = None
    user_location: UserLocation | None = None
    search_type: Literal["fast", "pro", "auto"] | None = None
    image_search_relevance_enhanced: bool | None = None


class MediaResponseOverrides(BaseModel):
    return_videos: bool | None = None
    return_images: bool | None = None


class MediaResponse(BaseModel):
    overrides: MediaResponseOverrides | None = None

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
