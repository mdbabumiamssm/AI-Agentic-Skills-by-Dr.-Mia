# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import re

import pytest
from langchain_core.utils import check_package_version


def test_check_package_version_pass() -> None:
    check_package_version("PyYAML", gte_version="5.4.1")


def test_check_package_version_fail() -> None:
    with pytest.raises(
        ValueError, match=re.escape("Expected PyYAML version to be < 5.4.1. Received ")
    ):
        check_package_version("PyYAML", lt_version="5.4.1")

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
