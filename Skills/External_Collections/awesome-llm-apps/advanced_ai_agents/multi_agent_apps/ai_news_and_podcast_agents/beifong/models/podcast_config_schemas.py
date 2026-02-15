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
from typing import Optional


class PodcastConfigBase(BaseModel):
    name: str
    prompt: str
    description: Optional[str] = None
    time_range_hours: int = Field(24, ge=1, le=168)
    limit_articles: int = Field(20, ge=5, le=50)
    is_active: bool = True
    tts_engine: str = "kokoro"
    language_code: str = "en"
    podcast_script_prompt: Optional[str] = None
    image_prompt: Optional[str] = None


class PodcastConfig(PodcastConfigBase):
    id: int
    created_at: Optional[str] = None
    updated_at: Optional[str] = None

    class Config:
        from_attributes = True


class PodcastConfigCreate(PodcastConfigBase):
    pass


class PodcastConfigUpdate(BaseModel):
    name: Optional[str] = None
    prompt: Optional[str] = None
    description: Optional[str] = None
    time_range_hours: Optional[int] = Field(None, ge=1, le=168)
    limit_articles: Optional[int] = Field(None, ge=5, le=50)
    is_active: Optional[bool] = None
    tts_engine: Optional[str] = None
    language_code: Optional[str] = None
    podcast_script_prompt: Optional[str] = None
    image_prompt: Optional[str] = None

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
