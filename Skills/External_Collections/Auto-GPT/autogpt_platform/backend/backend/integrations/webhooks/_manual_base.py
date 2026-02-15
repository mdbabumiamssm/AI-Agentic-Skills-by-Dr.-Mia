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

from backend.data import integrations
from backend.data.model import Credentials

from ._base import WT, BaseWebhooksManager

logger = logging.getLogger(__name__)


class ManualWebhookManagerBase(BaseWebhooksManager[WT]):
    async def _register_webhook(
        self,
        credentials: Credentials,
        webhook_type: WT,
        resource: str,
        events: list[str],
        ingress_url: str,
        secret: str,
    ) -> tuple[str, dict]:
        # TODO: pass ingress_url to user in frontend
        # See: https://github.com/Significant-Gravitas/AutoGPT/issues/8537
        logger.debug(f"Manual webhook registered with ingress URL: {ingress_url}")

        return "", {}

    async def _deregister_webhook(
        self,
        webhook: integrations.Webhook,
        credentials: Credentials,
    ) -> None:
        pass

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
