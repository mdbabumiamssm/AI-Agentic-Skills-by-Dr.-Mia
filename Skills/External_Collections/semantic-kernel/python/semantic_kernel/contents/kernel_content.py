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

from abc import ABC, abstractmethod
from typing import Annotated, Any, TypeVar

from pydantic import Field

from semantic_kernel.kernel_pydantic import KernelBaseModel

_T = TypeVar("_T", bound="KernelContent")


class KernelContent(KernelBaseModel, ABC):
    """Base class for all kernel contents."""

    # NOTE: if you wish to hold on to the inner content, you are responsible
    # for saving it before serializing the content/chat history as it won't be included.
    inner_content: Annotated[Any | None, Field(exclude=True)] = None
    ai_model_id: str | None = None
    metadata: dict[str, Any] = Field(default_factory=dict)

    @abstractmethod
    def __str__(self) -> str:
        """Return the string representation of the content."""
        pass

    @abstractmethod
    def to_element(self) -> Any:
        """Convert the instance to an Element."""
        pass

    @classmethod
    @abstractmethod
    def from_element(cls: type[_T], element: Any) -> _T:
        """Create an instance from an Element."""
        pass

    @abstractmethod
    def to_dict(self) -> dict[str, Any]:
        """Convert the instance to a dictionary."""
        pass

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
