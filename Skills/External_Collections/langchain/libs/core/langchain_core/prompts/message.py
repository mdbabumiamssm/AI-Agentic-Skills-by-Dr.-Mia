# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Message prompt templates."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Any

from langchain_core.load import Serializable
from langchain_core.utils.interactive_env import is_interactive_env

if TYPE_CHECKING:
    from langchain_core.messages import BaseMessage
    from langchain_core.prompts.chat import ChatPromptTemplate


class BaseMessagePromptTemplate(Serializable, ABC):
    """Base class for message prompt templates."""

    @classmethod
    def is_lc_serializable(cls) -> bool:
        """Return `True` as this class is serializable."""
        return True

    @classmethod
    def get_lc_namespace(cls) -> list[str]:
        """Get the namespace of the LangChain object.

        Returns:
            `["langchain", "prompts", "chat"]`
        """
        return ["langchain", "prompts", "chat"]

    @abstractmethod
    def format_messages(self, **kwargs: Any) -> list[BaseMessage]:
        """Format messages from kwargs. Should return a list of `BaseMessage` objects.

        Args:
            **kwargs: Keyword arguments to use for formatting.

        Returns:
            List of `BaseMessage` objects.
        """

    async def aformat_messages(self, **kwargs: Any) -> list[BaseMessage]:
        """Async format messages from kwargs.

        Args:
            **kwargs: Keyword arguments to use for formatting.

        Returns:
            List of `BaseMessage` objects.
        """
        return self.format_messages(**kwargs)

    @property
    @abstractmethod
    def input_variables(self) -> list[str]:
        """Input variables for this prompt template.

        Returns:
            List of input variables.
        """

    def pretty_repr(
        self,
        html: bool = False,  # noqa: FBT001,FBT002
    ) -> str:
        """Human-readable representation.

        Args:
            html: Whether to format as HTML.

        Returns:
            Human-readable representation.
        """
        raise NotImplementedError

    def pretty_print(self) -> None:
        """Print a human-readable representation."""
        print(self.pretty_repr(html=is_interactive_env()))  # noqa: T201

    def __add__(self, other: Any) -> ChatPromptTemplate:
        """Combine two prompt templates.

        Args:
            other: Another prompt template.

        Returns:
            Combined prompt template.
        """
        # Import locally to avoid circular import.
        from langchain_core.prompts.chat import ChatPromptTemplate  # noqa: PLC0415

        prompt = ChatPromptTemplate(messages=[self])
        return prompt.__add__(other)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
