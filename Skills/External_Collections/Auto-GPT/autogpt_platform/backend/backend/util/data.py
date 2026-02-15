# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import os
import pathlib
import sys


def get_frontend_path() -> pathlib.Path:
    if getattr(sys, "frozen", False):
        # The application is frozen
        datadir = pathlib.Path(os.path.dirname(sys.executable)) / "example_files"
    else:
        # The application is not frozen
        # Change this bit to match where you store your data files:
        filedir = os.path.dirname(__file__)
        datadir = pathlib.Path(filedir).parent.parent.parent / "example_files"
    return pathlib.Path(datadir)


def get_data_path() -> pathlib.Path:
    if getattr(sys, "frozen", False):
        # The application is frozen
        datadir = os.path.dirname(sys.executable)
    else:
        # The application is not frozen
        # Change this bit to match where you store your data files:
        filedir = os.path.dirname(__file__)
        datadir = pathlib.Path(filedir).parent.parent
    return pathlib.Path(datadir)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
