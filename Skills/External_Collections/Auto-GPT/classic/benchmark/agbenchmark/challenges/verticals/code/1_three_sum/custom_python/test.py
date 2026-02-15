# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

# pyright: reportMissingImports=false
from typing import List

from sample_code import three_sum


def test_three_sum(nums: List[int], target: int, expected_result: List[int]) -> None:
    result = three_sum(nums, target)
    print(result)
    assert (
        result == expected_result
    ), f"AssertionError: Expected the output to be {expected_result}"


if __name__ == "__main__":
    # test the trivial case with the first three numbers
    nums = [2, 7, 11, 15]
    target = 20
    expected_result = [0, 1, 2]
    test_three_sum(nums, target, expected_result)

    # test for ability to use zero and the same number twice
    nums = [2, 7, 0, 15, 12, 0]
    target = 2
    expected_result = [0, 2, 5]
    test_three_sum(nums, target, expected_result)

    # test for first and last index usage and negative numbers
    nums = [-6, 7, 11, 4]
    target = 9
    expected_result = [0, 2, 3]
    test_three_sum(nums, target, expected_result)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
