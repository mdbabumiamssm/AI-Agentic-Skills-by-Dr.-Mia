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
from abc import ABC, abstractmethod
from collections.abc import Awaitable, Callable
from typing import TypeVar

T = TypeVar("T")

logger: logging.Logger = logging.getLogger(__name__)


class RetryMechanismBase(ABC):
    """Base class for retry mechanisms."""

    @abstractmethod
    async def execute_with_retry(self, action: Callable[[], Awaitable[T]]) -> Awaitable[T]:
        """Executes the given action with retry logic.

        Args:
            action (Callable[[], Awaitable[T]]): The action to retry on exception.

        Returns:
            Awaitable[T]: An awaitable that will return the result of the action.
        """
        pass

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
