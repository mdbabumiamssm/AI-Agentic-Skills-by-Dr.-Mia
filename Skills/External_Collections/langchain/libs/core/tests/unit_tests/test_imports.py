# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import concurrent.futures
import importlib
import subprocess
from pathlib import Path


def test_importable_all() -> None:
    for path in Path("../core/langchain_core/").glob("*"):
        module_name = path.stem
        if (
            not module_name.startswith(".")
            and path.suffix != ".typed"
            and module_name != "pydantic_v1"
        ):
            module = importlib.import_module("langchain_core." + module_name)
            all_ = getattr(module, "__all__", [])
            for cls_ in all_:
                getattr(module, cls_)


def try_to_import(module_name: str) -> tuple[int, str]:
    """Try to import a module via subprocess."""
    module = importlib.import_module("langchain_core." + module_name)
    all_ = getattr(module, "__all__", [])
    for cls_ in all_:
        getattr(module, cls_)

    result = subprocess.run(
        ["python", "-c", f"import langchain_core.{module_name}"], check=True
    )
    return result.returncode, module_name


def test_importable_all_via_subprocess() -> None:
    """Test import in isolation.

    !!! note
        ImportErrors due to circular imports can be raised for one sequence of imports
        but not another.
    """
    module_names = []
    for path in Path("../core/langchain_core/").glob("*"):
        module_name = path.stem
        if (
            not module_name.startswith(".")
            and path.suffix != ".typed"
            and module_name != "pydantic_v1"
        ):
            module_names.append(module_name)

    with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
        futures = [
            executor.submit(try_to_import, module_name) for module_name in module_names
        ]
        for future in concurrent.futures.as_completed(futures):
            result = future.result()  # Will raise an exception if the callable raised
            code, module_name = result
            if code != 0:
                msg = f"Failed to import {module_name}."
                raise ValueError(msg)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
