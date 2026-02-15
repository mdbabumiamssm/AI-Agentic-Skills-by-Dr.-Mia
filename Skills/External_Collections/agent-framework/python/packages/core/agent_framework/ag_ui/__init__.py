# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

# Copyright (c) Microsoft. All rights reserved.

import importlib
from typing import Any

IMPORT_PATH = "agent_framework_ag_ui"
PACKAGE_NAME = "agent-framework-ag-ui"
_IMPORTS = [
    "__version__",
    "AgentFrameworkAgent",
    "add_agent_framework_fastapi_endpoint",
    "AGUIChatClient",
    "AGUIEventConverter",
    "AGUIHttpService",
    "ConfirmationStrategy",
    "DefaultConfirmationStrategy",
    "TaskPlannerConfirmationStrategy",
    "RecipeConfirmationStrategy",
    "DocumentWriterConfirmationStrategy",
]


def __getattr__(name: str) -> Any:
    if name in _IMPORTS:
        try:
            return getattr(importlib.import_module(IMPORT_PATH), name)
        except ModuleNotFoundError as exc:
            raise ModuleNotFoundError(
                f"The '{PACKAGE_NAME}' package is not installed, please do `pip install {PACKAGE_NAME}`"
            ) from exc
    raise AttributeError(f"Module {IMPORT_PATH} has no attribute {name}.")


def __dir__() -> list[str]:
    return _IMPORTS

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
