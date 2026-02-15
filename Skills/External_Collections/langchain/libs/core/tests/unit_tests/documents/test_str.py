# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from langchain_core.documents import Document


def test_str() -> None:
    assert str(Document(page_content="Hello, World!")) == "page_content='Hello, World!'"
    assert (
        str(Document(page_content="Hello, World!", metadata={"a": 3}))
        == "page_content='Hello, World!' metadata={'a': 3}"
    )


def test_repr() -> None:
    assert (
        repr(Document(page_content="Hello, World!"))
        == "Document(metadata={}, page_content='Hello, World!')"
    )
    assert (
        repr(Document(page_content="Hello, World!", metadata={"a": 3}))
        == "Document(metadata={'a': 3}, page_content='Hello, World!')"
    )

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
