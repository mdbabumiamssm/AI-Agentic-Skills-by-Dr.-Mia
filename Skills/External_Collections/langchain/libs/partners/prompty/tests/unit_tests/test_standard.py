# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import pytest
from pytest_benchmark.fixture import BenchmarkFixture  # type: ignore[import]

from langchain_prompty import create_chat_prompt


@pytest.mark.benchmark
def test_create_chat_prompt_init_time(benchmark: BenchmarkFixture) -> None:
    """Test create_chat_prompt initialization time."""

    def _create_chat_prompts() -> None:
        for _ in range(10):
            create_chat_prompt("Hello world")

    benchmark(_create_chat_prompts)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
