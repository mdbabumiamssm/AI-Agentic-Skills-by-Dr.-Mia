# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Tracers that record execution of LangChain runs."""

from typing import TYPE_CHECKING, Any

from langchain_core.tracers.langchain import LangChainTracer
from langchain_core.tracers.stdout import (
    ConsoleCallbackHandler,
    FunctionCallbackHandler,
)

from langchain_classic._api import create_importer
from langchain_classic.callbacks.tracers.logging import LoggingCallbackHandler

if TYPE_CHECKING:
    from langchain_community.callbacks.tracers.wandb import WandbTracer

# Create a way to dynamically look up deprecated imports.
# Used to consolidate logic for raising deprecation warnings and
# handling optional imports.
DEPRECATED_LOOKUP = {"WandbTracer": "langchain_community.callbacks.tracers.wandb"}

_import_attribute = create_importer(__file__, deprecated_lookups=DEPRECATED_LOOKUP)


def __getattr__(name: str) -> Any:
    """Look up attributes dynamically."""
    return _import_attribute(name)


__all__ = [
    "ConsoleCallbackHandler",
    "FunctionCallbackHandler",
    "LangChainTracer",
    "LoggingCallbackHandler",
    "WandbTracer",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
