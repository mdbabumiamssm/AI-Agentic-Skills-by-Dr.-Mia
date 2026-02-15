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


class VectorStoreException(KernelException):
    """Base class for all vector store exceptions."""

    pass


class VectorStoreInitializationException(VectorStoreException):
    """Class for all vector store initialization exceptions."""

    pass


class VectorStoreModelException(VectorStoreException):
    """Base class for all vector store model exceptions."""

    pass


class VectorStoreModelSerializationException(VectorStoreModelException):
    """An error occurred while serializing the vector store model."""

    pass


class VectorStoreModelDeserializationException(VectorStoreModelException):
    """An error occurred while deserializing the vector store model."""

    pass


class VectorStoreModelValidationError(VectorStoreModelException):
    """An error occurred while validating the vector store model."""

    pass


class VectorStoreMixinException(VectorStoreException):
    """Raised when a mixin is used without the VectorSearchBase Class."""

    pass


class VectorStoreOperationException(VectorStoreException):
    """An error occurred while performing an operation on the vector store record collection."""

    pass


class VectorStoreOperationNotSupportedException(VectorStoreOperationException):
    """An error occurred while performing an operation on the vector store record collection."""

    pass


class VectorSearchExecutionException(VectorStoreOperationException):
    """Raised when there is an error executing a VectorSearch function."""

    pass


class VectorSearchOptionsException(VectorStoreOperationException):
    """Raised when invalid options are given to a VectorSearch function."""

    pass


__all__ = [
    "VectorSearchExecutionException",
    "VectorSearchOptionsException",
    "VectorStoreException",
    "VectorStoreInitializationException",
    "VectorStoreMixinException",
    "VectorStoreModelDeserializationException",
    "VectorStoreModelException",
    "VectorStoreModelException",
    "VectorStoreModelSerializationException",
    "VectorStoreModelValidationError",
    "VectorStoreOperationException",
    "VectorStoreOperationNotSupportedException",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
