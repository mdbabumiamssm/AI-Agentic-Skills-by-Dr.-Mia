# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from typing import Type

from __module_name__.retrievers import __ModuleName__Retriever
from langchain_tests.integration_tests import (
    RetrieversIntegrationTests,
)


class Test__ModuleName__Retriever(RetrieversIntegrationTests):
    @property
    def retriever_constructor(self) -> Type[__ModuleName__Retriever]:
        """Get an empty vectorstore for unit tests."""
        return __ModuleName__Retriever

    @property
    def retriever_constructor_params(self) -> dict:
        return {"k": 2}

    @property
    def retriever_query_example(self) -> str:
        """Returns a str representing the "query" of an example retriever call."""
        return "example query"

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
