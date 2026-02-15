# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import os
from typing import Any

from langchain_core.utils import convert_to_secret_str
from perplexity import Perplexity


def initialize_client(values: dict[str, Any]) -> dict[str, Any]:
    """Initialize the Perplexity client."""
    pplx_api_key = (
        values.get("pplx_api_key")
        or os.environ.get("PPLX_API_KEY")
        or os.environ.get("PERPLEXITY_API_KEY")
        or ""
    )
    values["pplx_api_key"] = convert_to_secret_str(pplx_api_key)

    api_key = (
        values["pplx_api_key"].get_secret_value() if values["pplx_api_key"] else None
    )

    if not values.get("client"):
        values["client"] = Perplexity(api_key=api_key)

    return values

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
