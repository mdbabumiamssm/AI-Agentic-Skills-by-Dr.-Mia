# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from __future__ import annotations

from typing import Any, Optional

from pydantic import BaseModel, Field


class TaskRequestBody(BaseModel):
    input: str = Field(
        min_length=1,
        description="Input prompt for the task.",
        examples=["Write the words you receive to the file 'output.txt'."],
    )
    additional_input: Optional[dict[str, Any]] = Field(default_factory=dict)


class TaskEvalRequestBody(TaskRequestBody):
    eval_id: str

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
