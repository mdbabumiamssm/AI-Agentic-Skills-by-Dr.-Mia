# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from pathlib import Path

import pytest
from dotenv import load_dotenv

# Getting the absolute path of the current file's directory
ABS_PATH = Path(__file__).resolve().parent

# Getting the absolute path of the project's root directory
PROJECT_DIR = ABS_PATH.parent.parent


# Loading the .env file if it exists
def _load_env() -> None:
    dotenv_path = PROJECT_DIR / "tests" / "integration_tests" / ".env"
    if dotenv_path.exists():
        load_dotenv(dotenv_path)


_load_env()


@pytest.fixture(scope="module")
def test_dir() -> Path:
    return PROJECT_DIR / "tests" / "integration_tests"


# This fixture returns a string containing the path to the cassette directory for the
# current module
@pytest.fixture(scope="module")
def vcr_cassette_dir(request: pytest.FixtureRequest) -> str:
    module = Path(request.module.__file__)
    return str(module.parent / "cassettes" / module.stem)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
