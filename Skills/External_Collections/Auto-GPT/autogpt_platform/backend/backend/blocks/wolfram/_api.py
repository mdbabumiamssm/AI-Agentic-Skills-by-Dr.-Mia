# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from backend.sdk import APIKeyCredentials, Requests


async def llm_api_call(credentials: APIKeyCredentials, question: str) -> str:
    params = {"appid": credentials.api_key.get_secret_value(), "input": question}
    response = await Requests().get(
        "https://www.wolframalpha.com/api/v1/llm-api", params=params
    )
    if not response.ok:
        raise ValueError(f"API request failed: {response.status} {response.text()}")

    answer = response.text() if response.text() else ""

    return answer

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
