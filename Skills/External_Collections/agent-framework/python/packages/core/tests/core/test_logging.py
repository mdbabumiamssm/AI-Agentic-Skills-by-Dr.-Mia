# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

# Copyright (c) Microsoft. All rights reserved.


import pytest

from agent_framework import get_logger
from agent_framework.exceptions import AgentFrameworkException


def test_get_logger():
    """Test that the logger is created with the correct name."""
    logger = get_logger()
    assert logger.name == "agent_framework"


def test_get_logger_custom_name():
    """Test that the logger can be created with a custom name."""
    custom_name = "agent_framework.custom"
    logger = get_logger(custom_name)
    assert logger.name == custom_name


def test_get_logger_invalid_name():
    """Test that an exception is raised for an invalid logger name."""
    with pytest.raises(AgentFrameworkException):
        get_logger("invalid_name")


def test_log(caplog):
    """Test that the logger can log messages and adheres to the expected format."""
    logger = get_logger()
    with caplog.at_level("DEBUG"):
        logger.debug("This is a debug message")
        assert len(caplog.records) == 1
        record = caplog.records[0]
        assert record.levelname == "DEBUG"
        assert record.message == "This is a debug message"
        assert record.name == "agent_framework"
        assert record.pathname.endswith("test_logging.py")

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
