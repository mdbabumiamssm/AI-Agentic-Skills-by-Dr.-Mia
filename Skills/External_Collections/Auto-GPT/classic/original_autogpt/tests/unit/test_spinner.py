# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import time

from autogpt.app.spinner import Spinner

ALMOST_DONE_MESSAGE = "Almost done..."
PLEASE_WAIT = "Please wait..."


def test_spinner_initializes_with_default_values():
    """Tests that the spinner initializes with default values."""
    with Spinner() as spinner:
        assert spinner.message == "Loading..."
        assert spinner.delay == 0.1


def test_spinner_initializes_with_custom_values():
    """Tests that the spinner initializes with custom message and delay values."""
    with Spinner(message=PLEASE_WAIT, delay=0.2) as spinner:
        assert spinner.message == PLEASE_WAIT
        assert spinner.delay == 0.2


#
def test_spinner_stops_spinning():
    """Tests that the spinner starts spinning and stops spinning without errors."""
    with Spinner() as spinner:
        time.sleep(1)
    assert not spinner.running


def test_spinner_can_be_used_as_context_manager():
    """Tests that the spinner can be used as a context manager."""
    with Spinner() as spinner:
        assert spinner.running
    assert not spinner.running

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
