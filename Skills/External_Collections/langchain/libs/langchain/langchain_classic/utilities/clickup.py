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
    from langchain_community.utilities.clickup import (
        ClickupAPIWrapper,
        Component,
        CUList,
        Member,
        Space,
        Task,
        Team,
    )

# Create a way to dynamically look up deprecated imports.
# Used to consolidate logic for raising deprecation warnings and
# handling optional imports.
DEPRECATED_LOOKUP = {
    "Component": "langchain_community.utilities.clickup",
    "Task": "langchain_community.utilities.clickup",
    "CUList": "langchain_community.utilities.clickup",
    "Member": "langchain_community.utilities.clickup",
    "Team": "langchain_community.utilities.clickup",
    "Space": "langchain_community.utilities.clickup",
    "ClickupAPIWrapper": "langchain_community.utilities.clickup",
}

_import_attribute = create_importer(__package__, deprecated_lookups=DEPRECATED_LOOKUP)


def __getattr__(name: str) -> Any:
    """Look up attributes dynamically."""
    return _import_attribute(name)


__all__ = [
    "CUList",
    "ClickupAPIWrapper",
    "Component",
    "Member",
    "Space",
    "Task",
    "Team",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
