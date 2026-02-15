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

from typing import TypeVar

from samples.concepts.resources.utils import Colors, print_with_color
from semantic_kernel.data.vector import VectorSearchResult

_T = TypeVar("_T")


def print_record(result: VectorSearchResult[_T] | None = None, record: _T | None = None):
    if result:
        record = result.record
    print_with_color(f"  Found id: {record.id}", Colors.CGREEN)
    if result and result.score is not None:
        print_with_color(f"    Score: {result.score}", Colors.CWHITE)
    print_with_color(f"    Title: {record.title}", Colors.CWHITE)
    print_with_color(f"    Content: {record.content}", Colors.CWHITE)
    print_with_color(f"    Tag: {record.tag}", Colors.CWHITE)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
