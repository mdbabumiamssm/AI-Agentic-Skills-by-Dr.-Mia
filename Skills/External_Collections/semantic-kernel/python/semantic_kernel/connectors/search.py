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

import importlib

_IMPORTS = {
    "GoogleSearch": ".google",
    "GoogleSearchSettings": ".google",
    "GoogleSearchResult": ".google",
    "GoogleSearchResponse": ".google",
    "GoogleSearchInformation": ".google",
    "BraveSearch": ".brave",
    "BraveSettings": ".brave",
    "BraveWebPages": ".brave",
    "BraveWebPage": ".brave",
    "BraveSearchResponse": ".brave",
}


def __getattr__(name: str):
    if name in _IMPORTS:
        submod_name = _IMPORTS[name]
        module = importlib.import_module(submod_name, package=__name__)
        return getattr(module, name)
    raise AttributeError(f"module {__name__} has no attribute {name}")


def __dir__():
    return list(_IMPORTS.keys())

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
