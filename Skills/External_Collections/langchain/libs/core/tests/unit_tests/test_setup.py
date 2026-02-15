# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import time

import pytest
from blockbuster import BlockingError

from langchain_core import sys_info


async def test_blockbuster_setup() -> None:
    """Check if blockbuster is correctly setup."""
    # Blocking call outside of langchain_core is allowed.
    time.sleep(0.01)  # noqa: ASYNC251
    with pytest.raises(BlockingError):
        # Blocking call from langchain_core raises BlockingError.
        sys_info.print_sys_info()

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
