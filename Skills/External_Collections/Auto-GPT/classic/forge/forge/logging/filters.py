# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import logging


class BelowLevelFilter(logging.Filter):
    """Filter for logging levels below a certain threshold."""

    def __init__(self, below_level: int):
        super().__init__()
        self.below_level = below_level

    def filter(self, record: logging.LogRecord):
        return record.levelno < self.below_level

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
