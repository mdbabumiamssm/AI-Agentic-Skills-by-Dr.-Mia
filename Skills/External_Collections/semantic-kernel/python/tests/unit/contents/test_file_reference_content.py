# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

# Copyright (c) Microsoft. All rights reserved.

from xml.etree.ElementTree import Element

import pytest

from semantic_kernel.contents.file_reference_content import FileReferenceContent


def test_create_empty():
    file_reference = FileReferenceContent()
    assert file_reference.file_id is None


def test_create_file_id():
    file_reference = FileReferenceContent(file_id="12345")
    assert file_reference.file_id == "12345"


def test_update_file_id():
    file_reference = FileReferenceContent()
    file_reference.file_id = "12345"
    assert file_reference.file_id == "12345"


def test_to_str():
    file_reference = FileReferenceContent(file_id="12345")
    assert str(file_reference) == "FileReferenceContent(file_id=12345)"


def test_to_element():
    file_reference = FileReferenceContent(file_id="12345")
    element = file_reference.to_element()
    assert element.tag == "file_reference"
    assert element.get("file_id") == "12345"


def test_from_element():
    element = Element("FileReferenceContent")
    element.set("file_id", "12345")
    file_reference = FileReferenceContent.from_element(element)
    assert file_reference.file_id == "12345"


def test_to_dict_simple():
    file_reference = FileReferenceContent(file_id="12345")
    assert file_reference.to_dict() == {
        "file_id": "12345",
    }


@pytest.mark.parametrize(
    "file_reference",
    [
        pytest.param(FileReferenceContent(file_id="12345"), id="file_id"),
        pytest.param(FileReferenceContent(), id="empty"),
    ],
)
def test_element_roundtrip(file_reference):
    element = file_reference.to_element()
    new_file_reference = FileReferenceContent.from_element(element)
    assert new_file_reference == file_reference


@pytest.mark.parametrize(
    "file_reference",
    [
        pytest.param(FileReferenceContent(file_id="12345"), id="file_id"),
        pytest.param(FileReferenceContent(), id="empty"),
    ],
)
def test_to_dict(file_reference):
    expected_dict = {
        "file_id": file_reference.file_id,
    }
    assert file_reference.to_dict() == expected_dict

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
