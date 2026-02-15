# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from . import AutoGPTWolframAlphaSearch

plugin = AutoGPTWolframAlphaSearch()


def _wolframalpha_search(query: str) -> str | list[str]:
    res = ""
    try:
        ans = plugin.api.query(query)
        res = next(ans.results).text
    except Exception as e:
        return f"'_wolframalpha_search' on query: '{query}' raised exception: '{e}'"
    return res

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
