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
import warnings
from pathlib import Path

# Attempt to recursively import all modules in langchain
PKG_ROOT = Path(__file__).parent.parent.parent


def test_import_all() -> None:
    """Generate the public API for this package."""
    with warnings.catch_warnings():
        warnings.filterwarnings(action="ignore", category=UserWarning)
        library_code = PKG_ROOT / "langchain"
        for path in library_code.rglob("*.py"):
            # Calculate the relative path to the module
            module_name = path.relative_to(PKG_ROOT).with_suffix("").as_posix().replace("/", ".")
            if module_name.endswith("__init__"):
                # Without init
                module_name = module_name.rsplit(".", 1)[0]

            mod = importlib.import_module(module_name)

            all_attrs = getattr(mod, "__all__", [])

            for name in all_attrs:
                # Attempt to import the name from the module
                try:
                    obj = getattr(mod, name)
                    assert obj is not None
                except Exception as e:
                    msg = f"Could not import {module_name}.{name}"
                    raise AssertionError(msg) from e


def test_import_all_using_dir() -> None:
    """Generate the public API for this package."""
    library_code = PKG_ROOT / "langchain"
    for path in library_code.rglob("*.py"):
        # Calculate the relative path to the module
        module_name = path.relative_to(PKG_ROOT).with_suffix("").as_posix().replace("/", ".")
        if module_name.endswith("__init__"):
            # Without init
            module_name = module_name.rsplit(".", 1)[0]

        try:
            mod = importlib.import_module(module_name)
        except ModuleNotFoundError as e:
            msg = f"Could not import {module_name}"
            raise ModuleNotFoundError(msg) from e
        attributes = dir(mod)

        for name in attributes:
            if name.strip().startswith("_"):
                continue
            # Attempt to import the name from the module
            getattr(mod, name)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
