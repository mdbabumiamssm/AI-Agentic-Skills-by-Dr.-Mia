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

import importlib.metadata
from typing import Final

try:
    _version = importlib.metadata.version(__name__)
except importlib.metadata.PackageNotFoundError:
    _version = "0.0.0"  # Fallback for development mode
__version__: Final[str] = _version

from ._agents import *  # noqa: F403
from ._clients import *  # noqa: F403
from ._logging import *  # noqa: F403
from ._mcp import *  # noqa: F403
from ._memory import *  # noqa: F403
from ._middleware import *  # noqa: F403
from ._telemetry import *  # noqa: F403
from ._threads import *  # noqa: F403
from ._tools import *  # noqa: F403
from ._types import *  # noqa: F403
from ._workflows import *  # noqa: F403

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
