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

from semantic_kernel.connectors.ai import PromptExecutionSettings


def test_init():
    settings = PromptExecutionSettings()
    assert settings.service_id is None
    assert settings.extension_data == {}


def test_init_with_data():
    ext_data = {"test": "test"}
    settings = PromptExecutionSettings(service_id="test", extension_data=ext_data)
    assert settings.service_id == "test"
    assert settings.extension_data["test"] == "test"

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
