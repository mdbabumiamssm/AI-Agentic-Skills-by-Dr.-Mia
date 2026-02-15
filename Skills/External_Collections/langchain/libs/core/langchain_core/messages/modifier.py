# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Message responsible for deleting other messages."""

from typing import Any, Literal

from langchain_core.messages.base import BaseMessage


class RemoveMessage(BaseMessage):
    """Message responsible for deleting other messages."""

    type: Literal["remove"] = "remove"
    """The type of the message (used for serialization)."""

    def __init__(
        self,
        id: str,
        **kwargs: Any,
    ) -> None:
        """Create a RemoveMessage.

        Args:
            id: The ID of the message to remove.
            **kwargs: Additional fields to pass to the message.

        Raises:
            ValueError: If the 'content' field is passed in kwargs.

        """
        if kwargs.pop("content", None):
            msg = "RemoveMessage does not support 'content' field."
            raise ValueError(msg)

        super().__init__("", id=id, **kwargs)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
