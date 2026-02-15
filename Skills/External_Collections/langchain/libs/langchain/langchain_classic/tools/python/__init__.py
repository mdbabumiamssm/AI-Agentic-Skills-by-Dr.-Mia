# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from typing import Any


def __getattr__(_: str = "") -> Any:
    msg = (
        "This tool has been moved to langchain experiment. "
        "This tool has access to a python REPL. "
        "For best practices make sure to sandbox this tool. "
        "Read https://github.com/langchain-ai/langchain/blob/master/SECURITY.md "
        "To keep using this code as is, install langchain experimental and "
        "update relevant imports replacing 'langchain' with 'langchain_experimental'"
    )
    raise AttributeError(msg)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
