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
    from langchain_community.agent_toolkits.openapi.planner import (
        RequestsDeleteToolWithParsing,
        RequestsGetToolWithParsing,
        RequestsPatchToolWithParsing,
        RequestsPostToolWithParsing,
        RequestsPutToolWithParsing,
        create_openapi_agent,
    )

# Create a way to dynamically look up deprecated imports.
# Used to consolidate logic for raising deprecation warnings and
# handling optional imports.
DEPRECATED_LOOKUP = {
    "RequestsGetToolWithParsing": (
        "langchain_community.agent_toolkits.openapi.planner"
    ),
    "RequestsPostToolWithParsing": (
        "langchain_community.agent_toolkits.openapi.planner"
    ),
    "RequestsPatchToolWithParsing": (
        "langchain_community.agent_toolkits.openapi.planner"
    ),
    "RequestsPutToolWithParsing": (
        "langchain_community.agent_toolkits.openapi.planner"
    ),
    "RequestsDeleteToolWithParsing": (
        "langchain_community.agent_toolkits.openapi.planner"
    ),
    "create_openapi_agent": "langchain_community.agent_toolkits.openapi.planner",
}

_import_attribute = create_importer(__package__, deprecated_lookups=DEPRECATED_LOOKUP)


def __getattr__(name: str) -> Any:
    """Look up attributes dynamically."""
    return _import_attribute(name)


__all__ = [
    "RequestsDeleteToolWithParsing",
    "RequestsGetToolWithParsing",
    "RequestsPatchToolWithParsing",
    "RequestsPostToolWithParsing",
    "RequestsPutToolWithParsing",
    "create_openapi_agent",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
