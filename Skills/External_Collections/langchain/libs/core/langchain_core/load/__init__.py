# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""**Load** module helps with serialization and deserialization."""

from typing import TYPE_CHECKING

from langchain_core._import_utils import import_attr

if TYPE_CHECKING:
    from langchain_core.load.dump import dumpd, dumps
    from langchain_core.load.load import InitValidator, loads
    from langchain_core.load.serializable import Serializable

# Unfortunately, we have to eagerly import load from langchain_core/load/load.py
# eagerly to avoid a namespace conflict. We want users to still be able to use
# `from langchain_core.load import load` to get the load function, but
# the `from langchain_core.load.load import load` absolute import should also work.
from langchain_core.load.load import load

__all__ = (
    "InitValidator",
    "Serializable",
    "dumpd",
    "dumps",
    "load",
    "loads",
)

_dynamic_imports = {
    "dumpd": "dump",
    "dumps": "dump",
    "InitValidator": "load",
    "loads": "load",
    "Serializable": "serializable",
}


def __getattr__(attr_name: str) -> object:
    module_name = _dynamic_imports.get(attr_name)
    result = import_attr(attr_name, module_name, __spec__.parent)
    globals()[attr_name] = result
    return result


def __dir__() -> list[str]:
    return list(__all__)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
