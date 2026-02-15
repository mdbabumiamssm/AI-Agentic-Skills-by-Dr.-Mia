# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from typing import Optional

from pydantic import Field
from pydantic_settings import BaseSettings, SettingsConfigDict


class RateLimitSettings(BaseSettings):
    redis_host: str = Field(
        default="redis://localhost:6379",
        description="Redis host",
        validation_alias="REDIS_HOST",
    )

    redis_port: str = Field(
        default="6379", description="Redis port", validation_alias="REDIS_PORT"
    )

    redis_password: Optional[str] = Field(
        default=None,
        description="Redis password",
        validation_alias="REDIS_PASSWORD",
    )

    requests_per_minute: int = Field(
        default=60,
        description="Maximum number of requests allowed per minute per API key",
        validation_alias="RATE_LIMIT_REQUESTS_PER_MINUTE",
    )

    model_config = SettingsConfigDict(case_sensitive=True, extra="ignore")


RATE_LIMIT_SETTINGS = RateLimitSettings()

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
