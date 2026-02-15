# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from .scenex_plugin import SceneXplain

MOCK_API_KEY = "secret"
MOCK_IMAGE = "https://example.com/image.png"
MOCK_DESCRIPTION = "example description"


def test_describe_image(requests_mock):
    requests_mock.post(
        SceneXplain.API_ENDPOINT,
        json={
            "result": [
                {
                    "image": MOCK_IMAGE,
                    "text": MOCK_DESCRIPTION,
                }
            ]
        },
    )

    scenex = SceneXplain(MOCK_API_KEY)
    result = scenex.describe_image(
        image=MOCK_IMAGE,
        algorithm="Dune",
        features=[],
        languages=[],
    )

    # Check the results
    assert result == {
        "image": MOCK_IMAGE,
        "description": MOCK_DESCRIPTION,
    }

    # Check that the mocked functions were called with the correct arguments
    requests_mock.request_history[0].json() == {
        "data": [
            {
                "image": MOCK_IMAGE,
                "algorithm": "Dune",
                "features": [],
                "languages": [],
            }
        ]
    }

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
