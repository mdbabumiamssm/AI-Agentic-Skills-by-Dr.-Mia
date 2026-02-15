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


class ProcessException(KernelException):
    """Base class for all process exceptions."""

    pass


class ProcessInvalidConfigurationException(ProcessException):
    """An invalid configuration was provided for the process."""

    pass


class ProcessTargetFunctionNameMismatchException(ProcessException):
    """A message targeting a function has resulted in a different function becoming invocable."""

    pass


class ProcessFunctionNotFoundException(ProcessException):
    """A function was not found in the process."""

    pass


class ProcessEventUndefinedException(ProcessException):
    """An event was not defined in the process."""

    pass

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
