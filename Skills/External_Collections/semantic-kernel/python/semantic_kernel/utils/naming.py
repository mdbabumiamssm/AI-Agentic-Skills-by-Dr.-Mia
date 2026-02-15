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

import random
import string


def generate_random_ascii_name(length: int = 16) -> str:
    """Generate a series of random ASCII characters of the specified length.

    As example, plugin/function names can contain upper/lowercase letters, and underscores

    Args:
        length (int): The length of the string to generate.

    Returns:
        A string of random ASCII characters of the specified length.
    """
    letters = string.ascii_letters
    return "".join(random.choices(letters, k=length))  # nosec

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
