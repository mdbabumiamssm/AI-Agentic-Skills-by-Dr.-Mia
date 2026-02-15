# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from __future__ import annotations


class File:
    def __init__(self, name: str, content: list[str] | None = None) -> None:
        self.name = name
        self.content = "\n".join(content or [])

    def __eq__(self, __value: object, /) -> bool:
        if not isinstance(__value, File):
            return NotImplemented

        if self.name != __value.name:
            return False

        return self.content == __value.content

    def __hash__(self) -> int:
        return hash((self.name, self.content))

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
