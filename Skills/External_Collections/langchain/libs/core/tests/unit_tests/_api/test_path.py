# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from pathlib import Path

from langchain_core._api import path

HERE = Path(__file__).parent

ROOT = HERE.parent.parent.parent


def test_as_import_path() -> None:
    """Test that the path is converted to a LangChain import path."""
    # Verify that default paths are correct

    # if editable install, check directory structure
    if path.PACKAGE_DIR == ROOT / "langchain_core":
        assert path.PACKAGE_DIR == ROOT / "langchain_core"

    # Verify that as import path works correctly
    assert path.as_import_path(HERE, relative_to=ROOT) == "tests.unit_tests._api"
    assert (
        path.as_import_path(__file__, relative_to=ROOT)
        == "tests.unit_tests._api.test_path"
    )
    assert (
        path.as_import_path(__file__, suffix="create_agent", relative_to=ROOT)
        == "tests.unit_tests._api.test_path.create_agent"
    )

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
