# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import langchain_core
from langchain_core.callbacks.manager import _get_debug
from langchain_core.globals import get_debug, set_debug


def test_debug_is_settable_via_setter() -> None:
    previous_value = langchain_core.globals._debug
    previous_fn_reading = _get_debug()
    assert previous_value == previous_fn_reading

    # Flip the value of the flag.
    set_debug(not previous_value)

    new_value = langchain_core.globals._debug
    new_fn_reading = _get_debug()

    try:
        # We successfully changed the value of `debug`.
        assert new_value != previous_value

        # If we access `debug` via a function used elsewhere in langchain,
        # it also sees the same new value.
        assert new_value == new_fn_reading

        # If we access `debug` via `get_debug()` we also get the same value.
        assert new_value == get_debug()
    finally:
        # Make sure we don't alter global state, even if the test fails.
        # Always reset `debug` to the value it had before.
        set_debug(previous_value)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
