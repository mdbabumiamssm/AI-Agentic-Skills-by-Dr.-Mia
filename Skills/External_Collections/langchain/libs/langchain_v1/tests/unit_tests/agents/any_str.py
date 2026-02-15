# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import re


class AnyStr(str):
    __slots__ = ("prefix",)

    def __init__(self, prefix: str | re.Pattern = "") -> None:
        super().__init__()
        self.prefix = prefix

    def __eq__(self, other: object) -> bool:
        return isinstance(other, str) and (
            other.startswith(self.prefix)
            if isinstance(self.prefix, str)
            else self.prefix.match(other)
        )

    def __hash__(self) -> int:
        return hash((str(self), self.prefix))

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
