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


class SearchException(KernelException):
    """Base class for all Search related exceptions."""

    pass


class SearchResultEmptyError(SearchException):
    """Raised when there are no hits in the search results."""

    pass


class TextSearchException(SearchException):
    """An error occurred while executing a text search function."""

    pass


class TextSearchOptionsException(SearchException):
    """Raised when invalid options are given to a TextSearch function."""

    pass


__all__ = [
    "SearchException",
    "SearchResultEmptyError",
    "TextSearchException",
    "TextSearchOptionsException",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
