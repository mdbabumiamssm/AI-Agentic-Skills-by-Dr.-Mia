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
from semantic_kernel.exceptions.kernel_exceptions import KernelException


class FilterException(KernelException):
    """Base class for all filter exceptions."""

    pass


class FilterManagementException(FilterException):
    """An error occurred while adding or removing the filter to/from the kernel."""

    pass


__all__ = [
    "FilterException",
    "FilterManagementException",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
