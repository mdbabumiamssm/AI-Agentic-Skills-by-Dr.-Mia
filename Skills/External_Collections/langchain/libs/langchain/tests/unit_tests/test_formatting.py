# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Test formatting functionality."""

import pytest
from langchain_core.utils import formatter


def test_valid_formatting() -> None:
    """Test formatting works as expected."""
    template = "This is a {foo} test."
    output = formatter.format(template, foo="good")
    expected_output = "This is a good test."
    assert output == expected_output


def test_does_not_allow_args() -> None:
    """Test formatting raises error when args are provided."""
    template = "This is a {} test."
    with pytest.raises(
        ValueError,
        match="No arguments should be provided, "
        "everything should be passed as keyword arguments",
    ):
        formatter.format(template, "good")


def test_allows_extra_kwargs() -> None:
    """Test formatting allows extra keyword arguments."""
    template = "This is a {foo} test."
    output = formatter.format(template, foo="good", bar="oops")
    expected_output = "This is a good test."
    assert output == expected_output

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
