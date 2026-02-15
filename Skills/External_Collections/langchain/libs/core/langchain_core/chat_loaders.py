# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Chat loaders."""

from abc import ABC, abstractmethod
from collections.abc import Iterator

from langchain_core.chat_sessions import ChatSession


class BaseChatLoader(ABC):
    """Base class for chat loaders."""

    @abstractmethod
    def lazy_load(self) -> Iterator[ChatSession]:
        """Lazy load the chat sessions.

        Returns:
            An iterator of chat sessions.
        """

    def load(self) -> list[ChatSession]:
        """Eagerly load the chat sessions into memory.

        Returns:
            A list of chat sessions.
        """
        return list(self.lazy_load())

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
