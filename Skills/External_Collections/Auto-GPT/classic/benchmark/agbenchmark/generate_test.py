# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""
AGBenchmark's test discovery endpoint for Pytest.

This module is picked up by Pytest's *_test.py file matching pattern, and all challenge
classes in the module that conform to the `Test*` pattern are collected.
"""

import importlib
import logging
from itertools import chain

from agbenchmark.challenges.builtin import load_builtin_challenges
from agbenchmark.challenges.webarena import load_webarena_challenges

logger = logging.getLogger(__name__)

DATA_CATEGORY = {}

# Load challenges and attach them to this module
for challenge in chain(load_builtin_challenges(), load_webarena_challenges()):
    # Attach the Challenge class to this module so it can be discovered by pytest
    module = importlib.import_module(__name__)
    setattr(module, challenge.__name__, challenge)

    # Build a map of challenge names and their primary category
    DATA_CATEGORY[challenge.info.name] = challenge.info.category[0].value

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
