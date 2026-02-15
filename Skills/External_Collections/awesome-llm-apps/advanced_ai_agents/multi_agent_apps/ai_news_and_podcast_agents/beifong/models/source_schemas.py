# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from pydantic import BaseModel
from typing import Optional, List


class SourceFeed(BaseModel):
    id: int
    feed_url: str
    feed_type: str
    description: Optional[str] = None
    is_active: bool
    created_at: str
    last_crawled: Optional[str] = None


class SourceFeedCreate(BaseModel):
    feed_url: str
    feed_type: str = "main"
    description: Optional[str] = None
    is_active: bool = True


class SourceBase(BaseModel):
    name: str
    url: Optional[str] = None
    categories: Optional[List[str]] = []
    description: Optional[str] = None
    is_active: bool = True


class Source(SourceBase):
    id: int
    created_at: Optional[str] = None
    last_crawled: Optional[str] = None

    class Config:
        from_attributes = True


class SourceCreate(SourceBase):
    feeds: Optional[List[SourceFeedCreate]] = []


class SourceUpdate(BaseModel):
    name: Optional[str] = None
    url: Optional[str] = None
    categories: Optional[List[str]] = None
    description: Optional[str] = None
    is_active: Optional[bool] = None


class SourceWithFeeds(Source):
    feeds: List[SourceFeed] = []


class PaginatedSources(BaseModel):
    items: List[Source]
    total: int
    page: int
    per_page: int
    total_pages: int
    has_next: bool
    has_prev: bool


class Category(BaseModel):
    id: int
    name: str
    description: Optional[str] = None

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
