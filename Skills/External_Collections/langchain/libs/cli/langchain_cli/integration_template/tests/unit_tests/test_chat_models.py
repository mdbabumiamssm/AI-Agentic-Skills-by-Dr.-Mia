# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Test chat model integration."""

from typing import Type

from __module_name__.chat_models import Chat__ModuleName__
from langchain_tests.unit_tests import ChatModelUnitTests


class TestChat__ModuleName__Unit(ChatModelUnitTests):
    @property
    def chat_model_class(self) -> Type[Chat__ModuleName__]:
        return Chat__ModuleName__

    @property
    def chat_model_params(self) -> dict:
        # These should be parameters used to initialize your integration for testing
        return {
            "model": "bird-brain-001",
            "temperature": 0,
            "parrot_buffer_length": 50,
        }

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
