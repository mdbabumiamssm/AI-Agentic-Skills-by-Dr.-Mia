# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from dataclasses import dataclass

DEFAULT_USER_ID = "3e53486c-cf57-477e-ba2a-cb02dc828e1a"
DEFAULT_EMAIL = "default@example.com"


# Using dataclass here to avoid adding dependency on pydantic
@dataclass(frozen=True)
class User:
    user_id: str
    email: str
    phone_number: str
    role: str

    @classmethod
    def from_payload(cls, payload):
        return cls(
            user_id=payload["sub"],
            email=payload.get("email", ""),
            phone_number=payload.get("phone", ""),
            role=payload["role"],
        )

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
