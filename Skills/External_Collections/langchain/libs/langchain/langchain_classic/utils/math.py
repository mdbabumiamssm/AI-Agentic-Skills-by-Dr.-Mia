# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from typing import TYPE_CHECKING, Any

from langchain_classic._api import create_importer

if TYPE_CHECKING:
    from langchain_community.utils.math import (
        cosine_similarity,
        cosine_similarity_top_k,
    )

# Create a way to dynamically look up deprecated imports.
# Used to consolidate logic for raising deprecation warnings and
# handling optional imports.
# Not marked as deprecated since we may want to move the functionality
# into langchain as long as we're OK with numpy as the dependency.
_MODULE_LOOKUP = {
    "cosine_similarity": "langchain_community.utils.math",
    "cosine_similarity_top_k": "langchain_community.utils.math",
}

_import_attribute = create_importer(__package__, module_lookup=_MODULE_LOOKUP)


def __getattr__(name: str) -> Any:
    """Look up attributes dynamically."""
    return _import_attribute(name)


__all__ = [
    "cosine_similarity",
    "cosine_similarity_top_k",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
