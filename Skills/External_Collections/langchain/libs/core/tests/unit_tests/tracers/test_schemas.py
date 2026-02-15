# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from langchain_core.tracers import schemas
from langchain_core.tracers.schemas import __all__ as schemas_all


def test_public_api() -> None:
    """Test for changes in the public API."""
    expected_all = [
        "Run",
    ]

    assert sorted(schemas_all) == expected_all

    # Assert that the object is actually present in the schema module
    for module_name in expected_all:
        assert hasattr(schemas, module_name)
        assert getattr(schemas, module_name) is not None

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
