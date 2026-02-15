# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import logging

from fastapi import Request
from strenum import StrEnum

from backend.data import integrations
from backend.integrations.providers import ProviderName
from backend.sdk import Credentials

from ._manual_base import ManualWebhookManagerBase

logger = logging.getLogger(__name__)


class CompassWebhookType(StrEnum):
    TRANSCRIPTION = "transcription"
    TASK = "task"


class CompassWebhookManager(ManualWebhookManagerBase):
    PROVIDER_NAME = ProviderName.COMPASS
    WebhookType = CompassWebhookType

    @classmethod
    async def validate_payload(
        cls,
        webhook: integrations.Webhook,
        request: Request,
        credentials: Credentials | None,
    ) -> tuple[dict, str]:
        payload = await request.json()
        event_type = CompassWebhookType.TRANSCRIPTION  # currently the only type

        return payload, event_type

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
