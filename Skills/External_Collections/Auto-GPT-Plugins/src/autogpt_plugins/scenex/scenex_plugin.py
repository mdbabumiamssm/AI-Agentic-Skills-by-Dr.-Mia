# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from typing import List, Union

import requests

Algorithm = Union["Aqua", "Bolt", "Comet", "Dune", "Ember", "Flash"]


class SceneXplain:
    API_ENDPOINT = "https://us-central1-causal-diffusion.cloudfunctions.net/describe"

    def __init__(self, api_key):
        self._api_key = api_key

    def describe_image(
        self,
        image: str,
        algorithm: Algorithm = "Dune",
        features: List[str] = [],
        languages: List[str] = [],
    ) -> str:
        headers = {
            "x-api-key": f"token {self._api_key}",
            "content-type": "application/json",
        }

        payload = {
            "data": [
                {
                    "image": image,
                    "algorithm": algorithm,
                    "features": features,
                    "languages": languages,
                }
            ]
        }

        response = requests.post(self.API_ENDPOINT, headers=headers, json=payload)
        result = response.json().get("result", [])
        img = result[0] if result else {}

        return {"image": image, "description": img.get("text", "")}

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
