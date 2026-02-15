# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from .config import verify_settings
from .dependencies import (
    get_optional_user_id,
    get_user_id,
    requires_admin_user,
    requires_user,
)
from .helpers import add_auth_responses_to_openapi
from .models import User

__all__ = [
    "verify_settings",
    "get_user_id",
    "requires_admin_user",
    "requires_user",
    "get_optional_user_id",
    "add_auth_responses_to_openapi",
    "User",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
