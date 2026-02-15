# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import pytest

from .utils import remove_color_codes


@pytest.mark.parametrize(
    "raw_text, clean_text",
    [
        (
            "COMMAND = \x1b[36mbrowse_website\x1b[0m  "
            "ARGUMENTS = \x1b[36m{'url': 'https://www.google.com',"
            " 'question': 'What is the capital of France?'}\x1b[0m",
            "COMMAND = browse_website  "
            "ARGUMENTS = {'url': 'https://www.google.com',"
            " 'question': 'What is the capital of France?'}",
        ),
        (
            "{'Schaue dir meine Projekte auf github () an, als auch meine Webseiten': "
            "'https://github.com/Significant-Gravitas/AutoGPT,"
            " https://discord.gg/autogpt und https://twitter.com/Auto_GPT'}",
            "{'Schaue dir meine Projekte auf github () an, als auch meine Webseiten': "
            "'https://github.com/Significant-Gravitas/AutoGPT,"
            " https://discord.gg/autogpt und https://twitter.com/Auto_GPT'}",
        ),
        ("", ""),
        ("hello", "hello"),
        ("hello\x1b[31m world", "hello world"),
        ("\x1b[36mHello,\x1b[32m World!", "Hello, World!"),
        (
            "\x1b[1m\x1b[31mError:\x1b[0m\x1b[31m file not found",
            "Error: file not found",
        ),
    ],
)
def test_remove_color_codes(raw_text, clean_text):
    assert remove_color_codes(raw_text) == clean_text

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
