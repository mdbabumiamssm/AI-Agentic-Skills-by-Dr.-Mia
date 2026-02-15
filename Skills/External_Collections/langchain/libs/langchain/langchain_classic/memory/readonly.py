# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from typing import Any

from langchain_classic.base_memory import BaseMemory


class ReadOnlySharedMemory(BaseMemory):
    """Memory wrapper that is read-only and cannot be changed."""

    memory: BaseMemory

    @property
    def memory_variables(self) -> list[str]:
        """Return memory variables."""
        return self.memory.memory_variables

    def load_memory_variables(self, inputs: dict[str, Any]) -> dict[str, str]:
        """Load memory variables from memory."""
        return self.memory.load_memory_variables(inputs)

    def save_context(self, inputs: dict[str, Any], outputs: dict[str, str]) -> None:
        """Nothing should be saved or changed."""

    def clear(self) -> None:
        """Nothing to clear, got a memory like a vault."""

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
