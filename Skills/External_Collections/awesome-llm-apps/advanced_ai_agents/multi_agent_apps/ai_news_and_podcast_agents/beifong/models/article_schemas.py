# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from pydantic import BaseModel, ConfigDict
from typing import Optional, List, Dict, Any


class ArticleBase(BaseModel):
    title: str
    url: Optional[str] = None
    published_date: str
    summary: Optional[str] = None
    content: Optional[str] = None
    categories: Optional[List[str]] = []
    source_name: Optional[str] = None


class Article(ArticleBase):
    id: int
    metadata: Optional[Dict[str, Any]] = {}

    model_config = ConfigDict(from_attributes=True)


class PaginatedArticles(BaseModel):
    items: List[Article]
    total: int
    page: int
    per_page: int
    total_pages: int
    has_next: bool
    has_prev: bool
__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
