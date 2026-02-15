# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

# Copyright (c) Microsoft. All rights reserved.

import json


class MockResponse:
    def __init__(self, response, status=200):
        self._response = response
        self.status = status

    async def text(self):
        return self._response

    async def json(self):
        return self._response

    def raise_for_status(self):
        pass

    @property
    async def content(self):
        yield json.dumps(self._response).encode("utf-8")
        yield json.dumps({"done": True}).encode("utf-8")

    async def __aexit__(self, exc_type, exc, tb):
        pass

    async def __aenter__(self):
        return self

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
