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


class ContentException(KernelException):
    """Base class for all content exceptions."""

    pass


class ContentInitializationError(ContentException):
    """An error occurred while initializing the content."""

    pass


class ContentSerializationError(ContentException):
    """An error occurred while serializing the content."""

    pass


class ContentAdditionException(ContentException):
    """An error occurred while adding content."""

    pass


class FunctionCallInvalidNameException(ContentException):
    """An error occurred while validating the function name."""

    pass


class FunctionCallInvalidArgumentsException(ContentException):
    """An error occurred while validating the function arguments."""

    pass


class ChatHistoryReducerException(ContentException):
    """An error occurred while reducing chat history."""

    pass


__all__ = [
    "ChatHistoryReducerException",
    "ContentAdditionException",
    "ContentException",
    "ContentInitializationError",
    "ContentSerializationError",
    "FunctionCallInvalidArgumentsException",
    "FunctionCallInvalidNameException",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
