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
    from langchain_community.utilities.vertexai import (
        create_retry_decorator,
        get_client_info,
        init_vertexai,
        raise_vertex_import_error,
    )

# Create a way to dynamically look up deprecated imports.
# Used to consolidate logic for raising deprecation warnings and
# handling optional imports.
DEPRECATED_LOOKUP = {
    "create_retry_decorator": "langchain_community.utilities.vertexai",
    "raise_vertex_import_error": "langchain_community.utilities.vertexai",
    "init_vertexai": "langchain_community.utilities.vertexai",
    "get_client_info": "langchain_community.utilities.vertexai",
}

_import_attribute = create_importer(__package__, deprecated_lookups=DEPRECATED_LOOKUP)


def __getattr__(name: str) -> Any:
    """Look up attributes dynamically."""
    return _import_attribute(name)


__all__ = [
    "create_retry_decorator",
    "get_client_info",
    "init_vertexai",
    "raise_vertex_import_error",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
