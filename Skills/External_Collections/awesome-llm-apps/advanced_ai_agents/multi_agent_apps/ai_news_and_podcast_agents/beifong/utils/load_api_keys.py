# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import os
from pathlib import Path
from dotenv import load_dotenv


def load_api_key(key_name="OPENAI_API_KEY"):
    env_path = Path(".") / ".env"
    load_dotenv(dotenv_path=env_path)
    return os.environ.get(key_name)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
