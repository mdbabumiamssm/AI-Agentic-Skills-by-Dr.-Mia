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

from backend.sdk import Credentials, ManualWebhookManagerBase, Webhook

logger = logging.getLogger(__name__)


class GenericWebhookType(StrEnum):
    PLAIN = "plain"


class GenericWebhooksManager(ManualWebhookManagerBase):
    WebhookType = GenericWebhookType

    @classmethod
    async def validate_payload(
        cls, webhook: Webhook, request: Request, credentials: Credentials | None = None
    ) -> tuple[dict, str]:
        payload = await request.json()
        event_type = GenericWebhookType.PLAIN

        return payload, event_type

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
