# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

def __getattr__(_: str = "") -> None:
    """Raise an error on import since is deprecated."""
    msg = (
        "This module has been moved to langchain-experimental. "
        "For more details: https://github.com/langchain-ai/langchain/discussions/11352."
        "To access this code, install it with `pip install langchain-experimental`."
        "`from langchain_experimental.llm_bash.base "
        "import LLMBashChain`"
    )
    raise AttributeError(msg)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
