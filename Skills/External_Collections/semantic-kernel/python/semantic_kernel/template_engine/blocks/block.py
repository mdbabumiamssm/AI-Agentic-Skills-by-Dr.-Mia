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
from typing import ClassVar

from pydantic import field_validator

from semantic_kernel.kernel_pydantic import KernelBaseModel
from semantic_kernel.template_engine.blocks.block_types import BlockTypes

logger: logging.Logger = logging.getLogger(__name__)


class Block(KernelBaseModel):
    """A block."""

    type: ClassVar[BlockTypes] = BlockTypes.UNDEFINED
    content: str

    @field_validator("content", mode="before")
    @classmethod
    def content_strip(cls, content: str):
        """Strip the content of the block."""
        return content.strip()

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
