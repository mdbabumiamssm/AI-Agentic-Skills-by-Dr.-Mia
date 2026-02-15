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

from datetime import timedelta

from semantic_kernel.agents.open_ai.run_polling_options import RunPollingOptions


def test_get_polling_interval_below_threshold():
    options = RunPollingOptions()
    iteration_count = 1
    expected_interval = timedelta(milliseconds=250)
    assert options.get_polling_interval(iteration_count) == expected_interval


def test_get_polling_interval_at_threshold():
    options = RunPollingOptions()
    iteration_count = 2
    expected_interval = timedelta(milliseconds=250)
    assert options.get_polling_interval(iteration_count) == expected_interval


def test_get_polling_interval_above_threshold():
    options = RunPollingOptions()
    iteration_count = 3
    expected_interval = timedelta(seconds=1)
    assert options.get_polling_interval(iteration_count) == expected_interval


def test_get_polling_interval_custom_threshold():
    options = RunPollingOptions(run_polling_backoff_threshold=5)
    iteration_count = 4
    expected_interval = timedelta(milliseconds=250)
    assert options.get_polling_interval(iteration_count) == expected_interval

    iteration_count = 6
    expected_interval = timedelta(seconds=1)
    assert options.get_polling_interval(iteration_count) == expected_interval


def test_get_polling_interval_custom_intervals():
    options = RunPollingOptions(
        run_polling_interval=timedelta(milliseconds=500), run_polling_backoff=timedelta(seconds=2)
    )
    iteration_count = 1
    expected_interval = timedelta(milliseconds=500)
    assert options.get_polling_interval(iteration_count) == expected_interval

    iteration_count = 3
    expected_interval = timedelta(seconds=2)
    assert options.get_polling_interval(iteration_count) == expected_interval

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
