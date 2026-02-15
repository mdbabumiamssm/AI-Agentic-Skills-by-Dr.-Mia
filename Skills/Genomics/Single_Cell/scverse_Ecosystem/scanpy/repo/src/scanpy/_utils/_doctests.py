# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from __future__ import annotations

from collections.abc import Callable


def doctest_needs[F: Callable](mod: str) -> Callable[[F], F]:
    """Mark function with doctest dependency."""

    def decorator(func: F) -> F:
        func._doctest_needs = mod
        return func

    return decorator


def doctest_skip[F: Callable](reason: str) -> Callable[[F], F]:
    """Mark function so doctest is skipped."""
    if not reason:
        msg = "reason must not be empty"
        raise ValueError(msg)

    def decorator(func: F) -> F:
        func._doctest_skip_reason = reason
        return func

    return decorator


def doctest_internet[F: Callable](func: F) -> F:
    """Mark function so doctest gets the internet mark."""
    func._doctest_internet = True
    return func

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
