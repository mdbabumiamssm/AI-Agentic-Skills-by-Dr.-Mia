# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import inspect
from typing import cast


def is_caller_internal(depth: int = 2) -> bool:
    """Return whether the caller at `depth` of this function is internal."""
    try:
        frame = inspect.currentframe()
    except AttributeError:
        return False
    if frame is None:
        return False
    try:
        for _ in range(depth):
            frame = frame.f_back
            if frame is None:
                return False
        # Directly access the module name from the frame's global variables
        module_globals = frame.f_globals
        caller_module_name = cast("str", module_globals.get("__name__", ""))
        return caller_module_name.startswith("langchain")
    finally:
        del frame

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
