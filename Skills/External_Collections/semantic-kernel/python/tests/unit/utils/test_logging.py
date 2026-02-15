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

import logging

from semantic_kernel.utils.logging import setup_logging


def test_setup_logging():
    """Test that the logging is setup correctly."""
    setup_logging()

    root_logger = logging.getLogger()
    assert root_logger.handlers
    assert any(isinstance(handler, logging.StreamHandler) for handler in root_logger.handlers)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
