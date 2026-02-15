# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Test functionality related to prompt utils."""

from langchain_core.example_selectors import sorted_values


def test_sorted_vals() -> None:
    """Test sorted values from dictionary."""
    test_dict = {"key2": "val2", "key1": "val1"}
    expected_response = ["val1", "val2"]
    assert sorted_values(test_dict) == expected_response

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
