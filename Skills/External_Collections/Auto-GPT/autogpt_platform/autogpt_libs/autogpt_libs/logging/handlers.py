# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from __future__ import annotations

import json
import logging


class JsonFileHandler(logging.FileHandler):
    def format(self, record: logging.LogRecord) -> str:
        record.json_data = json.loads(record.getMessage())
        return json.dumps(getattr(record, "json_data"), ensure_ascii=False, indent=4)

    def emit(self, record: logging.LogRecord) -> None:
        with open(self.baseFilename, "w", encoding="utf-8") as f:
            f.write(self.format(record))

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
