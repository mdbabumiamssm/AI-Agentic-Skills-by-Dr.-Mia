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

from typing import Any

from langchain_classic._api import create_importer

# Code has been removed from the community package as well.
# We'll proxy to community package, which will raise an appropriate exception,
# but we'll not include this in __all__, so it won't be listed as importable.

_importer = create_importer(
    __package__,
    deprecated_lookups={"PythonREPL": "langchain_community.utilities.python"},
)


def __getattr__(name: str) -> Any:
    """Look up attributes dynamically."""
    return _importer(name)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
