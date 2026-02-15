# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from typing import Any, Dict, Optional

from pydantic import BaseModel


class Document(BaseModel):
    url: str
    relevance_score: float


class ApiResponse(BaseModel):
    answer: str
    documents: list[Document]
    success: bool


class GraphData(BaseModel):
    nodes: list[Dict[str, Any]]
    edges: list[Dict[str, Any]]
    graph_name: Optional[str] = None
    graph_description: Optional[str] = None


class Message(BaseModel):
    query: str
    response: str


class ChatRequest(BaseModel):
    query: str
    conversation_history: list[Message]
    message_id: str
    include_graph_data: bool = False
    graph_id: Optional[str] = None

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
