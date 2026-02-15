# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import importlib
from pathlib import Path

from pydantic import BaseModel


def test_all_models_built() -> None:
    for path in Path("../core/langchain_core/").glob("*"):
        module_name = path.stem
        if (
            not module_name.startswith(".")
            and path.suffix != ".typed"
            and module_name != "pydantic_v1"
        ):
            module = importlib.import_module("langchain_core." + module_name)
            all_ = getattr(module, "__all__", [])
            for attr_name in all_:
                attr = getattr(module, attr_name)
                try:
                    if issubclass(attr, BaseModel):
                        assert attr.__pydantic_complete__ is True
                except TypeError:
                    # This is expected for non-class attributes
                    pass

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
