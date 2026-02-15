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


class MemoryConnectorException(KernelException):
    """Base class for all memory connector exceptions."""

    pass


class MemoryConnectorConnectionException(MemoryConnectorException):
    """An error occurred while connecting to the memory connector."""

    pass


class MemoryConnectorInitializationError(MemoryConnectorException):
    """An error occurred while initializing the memory connector."""

    pass


class MemoryConnectorResourceNotFound(MemoryConnectorException):
    """The requested resource was not found in the memory connector."""

    pass


__all__ = [
    "MemoryConnectorConnectionException",
    "MemoryConnectorException",
    "MemoryConnectorInitializationError",
    "MemoryConnectorResourceNotFound",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
