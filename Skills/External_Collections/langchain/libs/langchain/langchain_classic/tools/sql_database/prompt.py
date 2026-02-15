# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""For backwards compatibility."""

from typing import TYPE_CHECKING, Any

from langchain_classic._api import create_importer

if TYPE_CHECKING:
    from langchain_community.tools.sql_database.prompt import QUERY_CHECKER


_importer = create_importer(
    __package__,
    deprecated_lookups={
        "QUERY_CHECKER": "langchain_community.tools.sql_database.prompt",
    },
)


def __getattr__(name: str) -> Any:
    """Look up attributes dynamically."""
    return _importer(name)


__all__ = ["QUERY_CHECKER"]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
