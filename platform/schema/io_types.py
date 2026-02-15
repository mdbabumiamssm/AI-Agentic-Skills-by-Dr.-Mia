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
from typing import List, Dict, Any, Optional
from datetime import datetime

class ToolCall(BaseModel):
    name: str
    arguments: Dict[str, Any]
    id: str

class LLMRequest(BaseModel):
    query: str
    system_instruction: Optional[str] = None
    tools: Optional[List[Dict[str, Any]]] = None
    temperature: float = 0.7
    max_tokens: int = 1024
    stop_sequences: Optional[List[str]] = None
    context: Optional[Dict[str, Any]] = {}

class LLMResponse(BaseModel):
    text: str
    tool_calls: Optional[List[ToolCall]] = None
    finish_reason: str = "stop" # stop, length, tool_calls, content_filter
    usage: Dict[str, int] = Field(default_factory=lambda: {"prompt_tokens": 0, "completion_tokens": 0})
    latency_ms: float = 0.0
    provider: str
    model: str
    timestamp: datetime = Field(default_factory=datetime.utcnow)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
