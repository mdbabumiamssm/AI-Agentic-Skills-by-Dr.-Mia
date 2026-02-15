# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""**Tracers** are classes for tracing runs."""

from typing import TYPE_CHECKING

from langchain_core._import_utils import import_attr

if TYPE_CHECKING:
    from langchain_core.tracers.base import BaseTracer
    from langchain_core.tracers.evaluation import EvaluatorCallbackHandler
    from langchain_core.tracers.langchain import LangChainTracer
    from langchain_core.tracers.log_stream import (
        LogStreamCallbackHandler,
        RunLog,
        RunLogPatch,
    )
    from langchain_core.tracers.schemas import Run
    from langchain_core.tracers.stdout import ConsoleCallbackHandler

__all__ = (
    "BaseTracer",
    "ConsoleCallbackHandler",
    "EvaluatorCallbackHandler",
    "LangChainTracer",
    "LogStreamCallbackHandler",
    "Run",
    "RunLog",
    "RunLogPatch",
)

_dynamic_imports = {
    "BaseTracer": "base",
    "EvaluatorCallbackHandler": "evaluation",
    "LangChainTracer": "langchain",
    "LogStreamCallbackHandler": "log_stream",
    "RunLog": "log_stream",
    "RunLogPatch": "log_stream",
    "Run": "schemas",
    "ConsoleCallbackHandler": "stdout",
}


def __getattr__(attr_name: str) -> object:
    module_name = _dynamic_imports.get(attr_name)
    result = import_attr(attr_name, module_name, __spec__.parent)
    globals()[attr_name] = result
    return result


def __dir__() -> list[str]:
    return list(__all__)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
