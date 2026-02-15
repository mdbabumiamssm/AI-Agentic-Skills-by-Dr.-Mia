# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""RunInfo class."""

from __future__ import annotations

from uuid import UUID

from pydantic import BaseModel


class RunInfo(BaseModel):
    """Class that contains metadata for a single execution of a Chain or model.

    Defined for backwards compatibility with older versions of langchain_core.

    This model will likely be deprecated in the future.

    Users can acquire the run_id information from callbacks or via run_id
    information present in the astream_event API (depending on the use case).
    """

    run_id: UUID
    """A unique identifier for the model or chain run."""

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
