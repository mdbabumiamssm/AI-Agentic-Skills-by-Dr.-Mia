# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from typing_extensions import override

from langchain_core.example_selectors import BaseExampleSelector


class DummyExampleSelector(BaseExampleSelector):
    def __init__(self) -> None:
        self.example: dict[str, str] | None = None

    def add_example(self, example: dict[str, str]) -> None:
        self.example = example

    @override
    def select_examples(self, input_variables: dict[str, str]) -> list[dict[str, str]]:
        return [input_variables]


async def test_aadd_example() -> None:
    selector = DummyExampleSelector()
    await selector.aadd_example({"foo": "bar"})
    assert selector.example == {"foo": "bar"}


async def test_aselect_examples() -> None:
    selector = DummyExampleSelector()
    examples = await selector.aselect_examples({"foo": "bar"})
    assert examples == [{"foo": "bar"}]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
