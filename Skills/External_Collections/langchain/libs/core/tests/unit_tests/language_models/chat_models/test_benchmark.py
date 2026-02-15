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
from itertools import cycle

from langchain_core.language_models import GenericFakeChatModel


def test_benchmark_model() -> None:
    """Add rate limiter."""
    tic = time.time()

    model = GenericFakeChatModel(
        messages=cycle(["hello", "world", "!"]),
    )

    for _ in range(1_000):
        model.invoke("foo")
    toc = time.time()

    # Verify that the time taken to run the loop is less than 1 seconds

    assert (toc - tic) < 1

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
