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
    from langchain_community.tools import APIOperation
    from langchain_community.tools.openapi.utils.api_models import (
        INVALID_LOCATION_TEMPL,
        PRIMITIVE_TYPES,
        SCHEMA_TYPE,
        SUPPORTED_LOCATIONS,
        APIProperty,
        APIPropertyBase,
        APIPropertyLocation,
        APIRequestBody,
        APIRequestBodyProperty,
    )

# Create a way to dynamically look up deprecated imports.
# Used to consolidate logic for raising deprecation warnings and
# handling optional imports.
DEPRECATED_LOOKUP = {
    "APIPropertyLocation": "langchain_community.tools.openapi.utils.api_models",
    "APIPropertyBase": "langchain_community.tools.openapi.utils.api_models",
    "APIProperty": "langchain_community.tools.openapi.utils.api_models",
    "APIRequestBodyProperty": "langchain_community.tools.openapi.utils.api_models",
    "APIRequestBody": "langchain_community.tools.openapi.utils.api_models",
    "APIOperation": "langchain_community.tools",
    "INVALID_LOCATION_TEMPL": "langchain_community.tools.openapi.utils.api_models",
    "SCHEMA_TYPE": "langchain_community.tools.openapi.utils.api_models",
    "PRIMITIVE_TYPES": "langchain_community.tools.openapi.utils.api_models",
    "SUPPORTED_LOCATIONS": "langchain_community.tools.openapi.utils.api_models",
}

_import_attribute = create_importer(__package__, deprecated_lookups=DEPRECATED_LOOKUP)


def __getattr__(name: str) -> Any:
    """Look up attributes dynamically."""
    return _import_attribute(name)


__all__ = [
    "INVALID_LOCATION_TEMPL",
    "PRIMITIVE_TYPES",
    "SCHEMA_TYPE",
    "SUPPORTED_LOCATIONS",
    "APIOperation",
    "APIProperty",
    "APIPropertyBase",
    "APIPropertyLocation",
    "APIRequestBody",
    "APIRequestBodyProperty",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
