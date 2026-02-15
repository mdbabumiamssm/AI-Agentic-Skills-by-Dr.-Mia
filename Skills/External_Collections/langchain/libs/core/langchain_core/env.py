# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Utilities for getting information about the runtime environment."""

import platform
from functools import lru_cache

from langchain_core import __version__


@lru_cache(maxsize=1)
def get_runtime_environment() -> dict:
    """Get information about the LangChain runtime environment.

    Returns:
        A dictionary with information about the runtime environment.
    """
    return {
        "library_version": __version__,
        "library": "langchain-core",
        "platform": platform.platform(),
        "runtime": "python",
        "runtime_version": platform.python_version(),
    }

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
