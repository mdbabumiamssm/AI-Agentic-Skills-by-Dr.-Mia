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

from langchain_core.language_models.base import _get_token_ids_default_method


class TestTokenCountingWithGPT2Tokenizer:
    def test_tokenization(self) -> None:
        # Check that the tokenization is consistent with the GPT-2 tokenizer
        assert _get_token_ids_default_method("This is a test") == [1212, 318, 257, 1332]

    def test_empty_token(self) -> None:
        assert len(_get_token_ids_default_method("")) == 0

    def test_multiple_tokens(self) -> None:
        assert len(_get_token_ids_default_method("a b c")) == 3

    def test_special_tokens(self) -> None:
        # test for consistency when the default tokenizer is changed
        assert len(_get_token_ids_default_method("a:b_c d")) == 6

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
