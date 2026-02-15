# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from typing import Any, Optional

from backend.util.request import Requests


class GetRequest:
    @classmethod
    async def get_request(
        cls, url: str, headers: Optional[dict] = None, json: bool = False
    ) -> Any:
        if headers is None:
            headers = {}
        response = await Requests().get(url, headers=headers)
        if json:
            return response.json()
        else:
            return response.text()

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
