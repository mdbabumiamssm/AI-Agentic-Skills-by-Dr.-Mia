# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Internal tracers used for stream_log and astream events implementations."""

import typing
from collections.abc import AsyncIterator, Iterator
from uuid import UUID

T = typing.TypeVar("T")


# THIS IS USED IN LANGGRAPH.
@typing.runtime_checkable
class _StreamingCallbackHandler(typing.Protocol[T]):
    """Types for streaming callback handlers.

    This is a common mixin that the callback handlers
    for both astream events and astream log inherit from.

    The `tap_output_aiter` method is invoked in some contexts
    to produce callbacks for intermediate results.
    """

    def tap_output_aiter(
        self, run_id: UUID, output: AsyncIterator[T]
    ) -> AsyncIterator[T]:
        """Used for internal astream_log and astream events implementations."""

    def tap_output_iter(self, run_id: UUID, output: Iterator[T]) -> Iterator[T]:
        """Used for internal astream_log and astream events implementations."""


__all__ = [
    "_StreamingCallbackHandler",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
