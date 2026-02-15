# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from langchain_core.utils import mustache
from pydantic import BaseModel

from .core import Invoker, Prompty, SimpleModel


class MustacheRenderer(Invoker):
    """Render a mustache template."""

    def __init__(self, prompty: Prompty) -> None:
        self.prompty = prompty

    def invoke(self, data: BaseModel) -> BaseModel:
        if not isinstance(data, SimpleModel):
            raise ValueError("Expected data to be an instance of SimpleModel")
        generated = mustache.render(self.prompty.content, data.item)
        return SimpleModel[str](item=generated)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
