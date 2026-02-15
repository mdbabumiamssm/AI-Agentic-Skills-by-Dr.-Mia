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

import click

logger = logging.getLogger(__name__)


def clean_input(prompt: str = ""):
    try:
        # ask for input, default when just pressing Enter is y
        logger.debug("Asking user via keyboard...")

        return click.prompt(
            text=prompt, prompt_suffix=" ", default="", show_default=False
        )
    except KeyboardInterrupt:
        logger.info("You interrupted AutoGPT")
        logger.info("Quitting...")
        exit(0)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
