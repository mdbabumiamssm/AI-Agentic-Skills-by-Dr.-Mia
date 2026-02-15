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

from forge.agent.forge_agent import ForgeAgent
from forge.agent_protocol.database.db import AgentDB
from forge.file_storage import FileStorageBackendName, get_storage

database_name = os.getenv("DATABASE_STRING")
workspace = get_storage(FileStorageBackendName.LOCAL, root_path=Path("workspace"))
database = AgentDB(database_name, debug_enabled=False)
agent = ForgeAgent(database=database, workspace=workspace)

app = agent.get_agent_app()

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
