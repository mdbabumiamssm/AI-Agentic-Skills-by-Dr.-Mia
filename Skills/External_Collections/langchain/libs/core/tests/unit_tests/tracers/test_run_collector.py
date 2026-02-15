# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Test the run collector."""

import uuid

from langchain_core.language_models import FakeListLLM
from langchain_core.tracers.context import collect_runs


def test_collect_runs() -> None:
    model = FakeListLLM(responses=["hello"])
    with collect_runs() as cb:
        model.invoke("hi")
        assert cb.traced_runs
        assert len(cb.traced_runs) == 1
        assert isinstance(cb.traced_runs[0].id, uuid.UUID)
        assert cb.traced_runs[0].inputs == {"prompts": ["hi"]}

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
