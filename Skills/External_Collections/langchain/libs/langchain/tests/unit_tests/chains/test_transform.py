# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Test transform chain."""

import re

import pytest

from langchain_classic.chains.transform import TransformChain


def dummy_transform(inputs: dict[str, str]) -> dict[str, str]:
    """Transform a dummy input for tests."""
    outputs = inputs
    outputs["greeting"] = f"{inputs['first_name']} {inputs['last_name']} says hello"
    del outputs["first_name"]
    del outputs["last_name"]
    return outputs


def test_transform_chain() -> None:
    """Test basic transform chain."""
    transform_chain = TransformChain(
        input_variables=["first_name", "last_name"],
        output_variables=["greeting"],
        transform=dummy_transform,
    )
    input_dict = {"first_name": "Leroy", "last_name": "Jenkins"}
    response = transform_chain(input_dict)
    expected_response = {"greeting": "Leroy Jenkins says hello"}
    assert response == expected_response


def test_transform_chain_bad_inputs() -> None:
    """Test basic transform chain."""
    transform_chain = TransformChain(
        input_variables=["first_name", "last_name"],
        output_variables=["greeting"],
        transform=dummy_transform,
    )
    input_dict = {"name": "Leroy", "last_name": "Jenkins"}
    with pytest.raises(
        ValueError, match=re.escape("Missing some input keys: {'first_name'}")
    ):
        _ = transform_chain(input_dict)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
