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

import logging
from abc import ABC

from pydantic import Field
from typing_extensions import deprecated

from semantic_kernel.kernel_pydantic import KernelBaseModel
from semantic_kernel.reliability.pass_through_without_retry import PassThroughWithoutRetry
from semantic_kernel.reliability.retry_mechanism_base import RetryMechanismBase

logger: logging.Logger = logging.getLogger(__name__)


class KernelReliabilityExtension(KernelBaseModel, ABC):
    """Kernel reliability extension."""

    retry_mechanism: RetryMechanismBase = Field(
        default_factory=PassThroughWithoutRetry,
        exclude=True,
        deprecated=deprecated("retry_mechanism is deprecated; This property doesn't have any effect on the kernel."),
    )

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
