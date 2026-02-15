# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from zerobouncesdk import ZBValidateResponse, ZeroBounce


class ZeroBounceClient:
    def __init__(self, api_key: str):
        self.api_key = api_key
        self.client = ZeroBounce(api_key)

    def validate_email(self, email: str, ip_address: str) -> ZBValidateResponse:
        return self.client.validate(email, ip_address)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
