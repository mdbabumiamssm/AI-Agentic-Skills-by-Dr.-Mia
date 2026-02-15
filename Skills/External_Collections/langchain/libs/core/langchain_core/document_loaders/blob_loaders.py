# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Schema for Blobs and Blob Loaders.

The goal is to facilitate decoupling of content loading from content parsing code.

In addition, content loading code should provide a lazy loading interface by default.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

# Re-export Blob and PathLike for backwards compatibility
from langchain_core.documents.base import Blob, PathLike

if TYPE_CHECKING:
    from collections.abc import Iterable


class BlobLoader(ABC):
    """Abstract interface for blob loaders implementation.

    Implementer should be able to load raw content from a storage system according
    to some criteria and return the raw content lazily as a stream of blobs.
    """

    @abstractmethod
    def yield_blobs(
        self,
    ) -> Iterable[Blob]:
        """A lazy loader for raw data represented by LangChain's `Blob` object.

        Returns:
            A generator over blobs
        """


# Re-export Blob and Pathlike for backwards compatibility
__all__ = ["Blob", "BlobLoader", "PathLike"]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
