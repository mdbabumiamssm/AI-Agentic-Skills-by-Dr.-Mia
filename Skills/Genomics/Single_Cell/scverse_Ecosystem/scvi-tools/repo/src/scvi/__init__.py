# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""scvi-tools."""

# Set default logging handler to avoid logging with logging.lastResort logger.
import logging

from ._constants import REGISTRY_KEYS
from ._settings import settings

# this import needs to come after prior imports to prevent circular import
from . import data, model, external, utils, criticism

from importlib.metadata import version

package_name = "scvi-tools"
__version__ = version(package_name)

settings.verbosity = logging.INFO

# Jax sets the root logger, this prevents double output.
scvi_logger = logging.getLogger("scvi")
scvi_logger.propagate = False


__all__ = [
    "settings",
    "REGISTRY_KEYS",
    "data",
    "model",
    "external",
    "utils",
    "criticism",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
