# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from typing import Any, Dict, List, Optional

from pydantic import BaseModel, Field


class AutoModRequest(BaseModel):
    """Request model for AutoMod API"""

    type: str = Field(..., description="Content type - 'text', 'image', 'video'")
    content: str = Field(..., description="The content to moderate")
    metadata: Optional[Dict[str, Any]] = Field(
        default=None, description="Additional context about the content"
    )


class ModerationResult(BaseModel):
    """Individual moderation result"""

    decision: str = Field(
        ..., description="Moderation decision: 'approved', 'rejected', 'flagged'"
    )
    reason: Optional[str] = Field(default=None, description="Reason for the decision")


class AutoModResponse(BaseModel):
    """Response model for AutoMod API"""

    success: bool = Field(..., description="Whether the request was successful")
    content_id: str = Field(
        ..., description="Unique reference ID for this moderation request"
    )
    status: str = Field(
        ..., description="Overall status: 'approved', 'rejected', 'flagged', 'pending'"
    )
    moderation_results: List[ModerationResult] = Field(
        default_factory=list, description="List of moderation results"
    )


class ModerationConfig(BaseModel):
    """Configuration for AutoMod integration"""

    enabled: bool = Field(default=True, description="Whether moderation is enabled")
    api_url: str = Field(default="", description="AutoMod API base URL")
    api_key: str = Field(..., description="AutoMod API key")
    timeout: int = Field(default=30, description="Request timeout in seconds")
    retry_attempts: int = Field(default=3, description="Number of retry attempts")
    retry_delay: float = Field(
        default=1.0, description="Delay between retries in seconds"
    )
    fail_open: bool = Field(
        default=False,
        description="If True, allow execution to continue if moderation fails",
    )
    moderate_inputs: bool = Field(
        default=True, description="Whether to moderate block inputs"
    )
    moderate_outputs: bool = Field(
        default=True, description="Whether to moderate block outputs"
    )

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
