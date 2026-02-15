# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import logging
import sys
import uuid

import pytest

from langchain_classic.callbacks.tracers import LoggingCallbackHandler


def test_logging(
    caplog: pytest.LogCaptureFixture,
    capsys: pytest.CaptureFixture[str],
) -> None:
    # Set up a Logger and a handler so we can check the Logger's handlers work too
    logger = logging.getLogger("test_logging")
    logger.setLevel(logging.INFO)
    logger.addHandler(logging.StreamHandler(sys.stdout))

    handler = LoggingCallbackHandler(logger, extra={"test": "test_extra"})
    handler.on_text("test", run_id=uuid.uuid4())

    # Assert logging actually took place
    assert len(caplog.record_tuples) == 1
    record = caplog.records[0]
    assert record.name == logger.name
    assert record.levelno == logging.INFO
    assert (
        record.msg == "\x1b[36;1m\x1b[1;3m[text]\x1b[0m \x1b[1mNew text:\x1b[0m\ntest"
    )
    # Check the extra shows up
    assert record.test == "test_extra"  # type: ignore[attr-defined]

    # Assert log handlers worked
    cap_result = capsys.readouterr()
    assert (
        cap_result.out
        == "\x1b[36;1m\x1b[1;3m[text]\x1b[0m \x1b[1mNew text:\x1b[0m\ntest\n"
    )

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
