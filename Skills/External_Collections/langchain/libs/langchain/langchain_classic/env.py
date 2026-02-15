# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import platform
from functools import lru_cache


@lru_cache(maxsize=1)
def get_runtime_environment() -> dict:
    """Get information about the LangChain runtime environment."""
    # Lazy import to avoid circular imports
    from langchain_classic import __version__

    return {
        "library_version": __version__,
        "library": "langchain-classic",
        "platform": platform.platform(),
        "runtime": "python",
        "runtime_version": platform.python_version(),
    }

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
