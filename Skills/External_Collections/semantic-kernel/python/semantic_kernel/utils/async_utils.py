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

import asyncio
from collections.abc import Callable
from functools import partial
from typing import Any


async def run_in_executor(executor: Any, func: Callable, *args, **kwargs) -> Any:
    """Run a function in an executor."""
    return await asyncio.get_event_loop().run_in_executor(executor, partial(func, *args, **kwargs))

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
