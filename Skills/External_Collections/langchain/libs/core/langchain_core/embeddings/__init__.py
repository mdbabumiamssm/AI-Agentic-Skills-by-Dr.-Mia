# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Embeddings."""

from typing import TYPE_CHECKING

from langchain_core._import_utils import import_attr

if TYPE_CHECKING:
    from langchain_core.embeddings.embeddings import Embeddings
    from langchain_core.embeddings.fake import (
        DeterministicFakeEmbedding,
        FakeEmbeddings,
    )

__all__ = ("DeterministicFakeEmbedding", "Embeddings", "FakeEmbeddings")

_dynamic_imports = {
    "Embeddings": "embeddings",
    "DeterministicFakeEmbedding": "fake",
    "FakeEmbeddings": "fake",
}


def __getattr__(attr_name: str) -> object:
    module_name = _dynamic_imports.get(attr_name)
    result = import_attr(attr_name, module_name, __spec__.parent)
    globals()[attr_name] = result
    return result


def __dir__() -> list[str]:
    return list(__all__)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
