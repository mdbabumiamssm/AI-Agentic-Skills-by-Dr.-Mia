# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Tests for the InMemoryStore class."""

import pytest
from langchain_core.stores import InMemoryStore
from typing_extensions import override

from langchain_tests.integration_tests.base_store import (
    BaseStoreAsyncTests,
    BaseStoreSyncTests,
)


class TestInMemoryStore(BaseStoreSyncTests[str]):
    @pytest.fixture
    @override
    def three_values(self) -> tuple[str, str, str]:
        return "foo", "bar", "buzz"

    @pytest.fixture
    @override
    def kv_store(self) -> InMemoryStore:
        return InMemoryStore()


class TestInMemoryStoreAsync(BaseStoreAsyncTests[str]):
    @pytest.fixture
    @override
    def three_values(self) -> tuple[str, str, str]:
        return "foo", "bar", "buzz"

    @pytest.fixture
    @override
    async def kv_store(self) -> InMemoryStore:
        return InMemoryStore()

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
