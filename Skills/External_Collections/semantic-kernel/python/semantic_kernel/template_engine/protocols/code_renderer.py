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

from abc import abstractmethod
from typing import TYPE_CHECKING, Protocol, runtime_checkable

if TYPE_CHECKING:
    from semantic_kernel import Kernel
    from semantic_kernel.functions.kernel_arguments import KernelArguments


@runtime_checkable
class CodeRenderer(Protocol):
    """Protocol for dynamic code blocks that need async IO to be rendered."""

    @abstractmethod
    async def render_code(self, kernel: "Kernel", arguments: "KernelArguments") -> str:
        """Render the block using the given context.

        :param context: kernel execution context
        :return: Rendered content
        """

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
