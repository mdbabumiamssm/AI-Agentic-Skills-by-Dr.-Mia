# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Helper functions for managing the LangChain API.

This module is only relevant for LangChain developers, not for users.

!!! warning

    This module and its submodules are for internal use only. Do not use them in your
    own code.  We may change the API at any time with no warning.

"""

from .deprecation import (
    LangChainDeprecationWarning,
    deprecated,
    suppress_langchain_deprecation_warning,
    surface_langchain_deprecation_warnings,
    warn_deprecated,
)
from .module_import import create_importer

__all__ = [
    "LangChainDeprecationWarning",
    "create_importer",
    "deprecated",
    "suppress_langchain_deprecation_warning",
    "surface_langchain_deprecation_warnings",
    "warn_deprecated",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
