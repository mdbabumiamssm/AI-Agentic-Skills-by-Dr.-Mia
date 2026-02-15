# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from unittest.mock import patch

import pytest
from forge.config.ai_directives import AIDirectives
from forge.config.ai_profile import AIProfile

from autogpt.app.config import AppConfig
from autogpt.app.setup import (
    apply_overrides_to_ai_settings,
    interactively_revise_ai_settings,
)


@pytest.mark.asyncio
async def test_apply_overrides_to_ai_settings():
    ai_profile = AIProfile(ai_name="Test AI", ai_role="Test Role")
    directives = AIDirectives(
        resources=["Resource1"],
        constraints=["Constraint1"],
        best_practices=["BestPractice1"],
    )

    apply_overrides_to_ai_settings(
        ai_profile,
        directives,
        override_name="New AI",
        override_role="New Role",
        replace_directives=True,
        resources=["NewResource"],
        constraints=["NewConstraint"],
        best_practices=["NewBestPractice"],
    )

    assert ai_profile.ai_name == "New AI"
    assert ai_profile.ai_role == "New Role"
    assert directives.resources == ["NewResource"]
    assert directives.constraints == ["NewConstraint"]
    assert directives.best_practices == ["NewBestPractice"]


@pytest.mark.asyncio
async def test_interactively_revise_ai_settings(config: AppConfig):
    ai_profile = AIProfile(ai_name="Test AI", ai_role="Test Role")
    directives = AIDirectives(
        resources=["Resource1"],
        constraints=["Constraint1"],
        best_practices=["BestPractice1"],
    )

    user_inputs = [
        "n",
        "New AI",
        "New Role",
        "NewConstraint",
        "",
        "NewResource",
        "",
        "NewBestPractice",
        "",
        "y",
    ]
    with patch("autogpt.app.setup.clean_input", side_effect=user_inputs):
        ai_profile, directives = await interactively_revise_ai_settings(
            ai_profile, directives, config
        )

    assert ai_profile.ai_name == "New AI"
    assert ai_profile.ai_role == "New Role"
    assert directives.resources == ["NewResource"]
    assert directives.constraints == ["NewConstraint"]
    assert directives.best_practices == ["NewBestPractice"]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
