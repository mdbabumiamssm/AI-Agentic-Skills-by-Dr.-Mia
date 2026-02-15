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

from collections.abc import Sequence
from typing import TypeVar, Union

from semantic_kernel.services.ai_service_client_base import AIServiceClientBase

AI_SERVICE_CLIENT_TYPE = TypeVar("AI_SERVICE_CLIENT_TYPE", bound=AIServiceClientBase)

T = TypeVar("T")

OneOrMany = Union[T, Sequence[T]]
OneOrList = Union[T, list[T]]
OptionalOneOrMany = Union[T, Sequence[T], None]
OptionalOneOrList = Union[T, list[T], None]

__all__ = ["AI_SERVICE_CLIENT_TYPE", "OneOrList", "OneOrMany", "OptionalOneOrList", "OptionalOneOrMany"]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
