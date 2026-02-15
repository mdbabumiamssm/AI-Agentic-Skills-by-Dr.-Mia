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
from collections.abc import Callable
from functools import wraps


def error_on_missing_dependencies(*modules):
    missing_modules = []
    for module in modules:
        try:
            importlib.import_module(module)
        except ImportError:
            missing_modules.append(module)
    if len(missing_modules) > 0:
        raise ModuleNotFoundError(f"Please install {missing_modules} to use this functionality.")


def dependencies(*modules) -> Callable:
    """Decorator to check for dependencies."""

    def decorator(fn: Callable) -> Callable:
        @wraps(fn)
        def wrapper(*args, **kwargs):
            error_on_missing_dependencies(*modules)
            return fn(*args, **kwargs)

        return wrapper

    return decorator


def is_package_installed(module):
    try:
        importlib.import_module(module)  # Try importing the package
        return True  # If no exception is raised, the package is installed
    except ImportError:
        return False  # If ImportError is raised, the package is not installed

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
