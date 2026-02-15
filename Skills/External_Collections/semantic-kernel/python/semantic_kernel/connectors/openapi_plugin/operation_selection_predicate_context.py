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


class OperationSelectionPredicateContext:
    """The context for the operation selection predicate."""

    def __init__(self, operation_id: str, path: str, method: str, description: str | None = None):
        """Initialize the operation selection predicate context."""
        self.operation_id = operation_id
        self.path = path
        self.method = method
        self.description = description

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
