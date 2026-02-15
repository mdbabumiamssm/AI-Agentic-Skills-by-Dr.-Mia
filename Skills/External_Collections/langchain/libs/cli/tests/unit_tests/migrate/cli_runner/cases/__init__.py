# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from tests.unit_tests.migrate.cli_runner.case import Case
from tests.unit_tests.migrate.cli_runner.cases import imports
from tests.unit_tests.migrate.cli_runner.file import File
from tests.unit_tests.migrate.cli_runner.folder import Folder

cases = [
    Case(
        name="empty",
        source=File("__init__.py", content=[]),
        expected=File("__init__.py", content=[]),
    ),
    *imports.cases,
]
before = Folder("project", *[case.source for case in cases])
expected = Folder("project", *[case.expected for case in cases])

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
