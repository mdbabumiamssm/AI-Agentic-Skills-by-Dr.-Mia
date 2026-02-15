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

from typing import ClassVar

from pydantic import SecretStr

from semantic_kernel.kernel_pydantic import KernelBaseSettings
from semantic_kernel.utils.feature_stage_decorator import experimental


@experimental
class AstraDBSettings(KernelBaseSettings):
    """AstraDB model settings.

    Settings for AstraDB connection:
    - app_token: SecretStr | None - AstraDB token (Env var ASTRADB_APP_TOKEN)
    - db_id: str | None - AstraDB database ID (Env var ASTRADB_DB_ID)
    - region: str | None - AstraDB region (Env var ASTRADB_REGION)
    - keyspace: str | None - AstraDB keyspace (Env var ASTRADB_KEYSPACE)
    """

    env_prefix: ClassVar[str] = "ASTRADB_"

    app_token: SecretStr
    db_id: str
    region: str
    keyspace: str

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
