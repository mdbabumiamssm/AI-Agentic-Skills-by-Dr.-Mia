# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

# Import the provider builder to ensure it's registered
from backend.sdk.registry import AutoRegistry

from .triggers import GenericWebhookTriggerBlock, generic_webhook

# Ensure the SDK registry is patched to include our webhook manager
AutoRegistry.patch_integrations()

__all__ = ["GenericWebhookTriggerBlock", "generic_webhook"]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
