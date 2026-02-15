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
from collections.abc import AsyncGenerator, AsyncIterable, Sequence
from typing import TypeVar

_T = TypeVar("_T")


async def desync_list(sync_list: Sequence[_T]) -> AsyncIterable[_T]:  # noqa: RUF029
    """De synchronize a list of synchronous objects."""
    for x in sync_list:
        yield x


async def empty_generator() -> AsyncGenerator[_T, None]:
    """An empty generator, can be used to return an empty generator."""
    if False:
        yield None
    await asyncio.sleep(0)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
