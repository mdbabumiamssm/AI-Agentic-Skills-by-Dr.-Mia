# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""GitHub utilities."""

from __future__ import annotations

import http.client
import json


def list_packages(*, contains: str | None = None) -> list[str]:
    """List all packages in the langchain repository templates directory.

    Args:
        contains: Optional substring that the package name must contain.

    Returns:
        A list of package names.
    """
    conn = http.client.HTTPSConnection("api.github.com")
    try:
        headers = {
            "Accept": "application/vnd.github+json",
            "X-GitHub-Api-Version": "2022-11-28",
            "User-Agent": "langchain-cli",
        }

        conn.request(
            "GET",
            "/repos/langchain-ai/langchain/contents/templates",
            headers=headers,
        )
        res = conn.getresponse()

        res_str = res.read()

        data = json.loads(res_str)
        package_names = [
            p["name"] for p in data if p["type"] == "dir" and p["name"] != "docs"
        ]
        return (
            [p for p in package_names if contains in p] if contains else package_names
        )
    finally:
        conn.close()

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
