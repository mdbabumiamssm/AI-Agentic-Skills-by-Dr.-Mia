# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Vector stores."""

from typing import TYPE_CHECKING

from langchain_core._import_utils import import_attr

if TYPE_CHECKING:
    from langchain_core.vectorstores.base import VST, VectorStore, VectorStoreRetriever
    from langchain_core.vectorstores.in_memory import InMemoryVectorStore

__all__ = (
    "VST",
    "InMemoryVectorStore",
    "VectorStore",
    "VectorStoreRetriever",
)

_dynamic_imports = {
    "VectorStore": "base",
    "VST": "base",
    "VectorStoreRetriever": "base",
    "InMemoryVectorStore": "in_memory",
}


def __getattr__(attr_name: str) -> object:
    module_name = _dynamic_imports.get(attr_name)
    result = import_attr(attr_name, module_name, __spec__.parent)
    globals()[attr_name] = result
    return result


def __dir__() -> list[str]:
    return list(__all__)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
