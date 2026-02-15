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

from semantic_kernel import Kernel
from semantic_kernel.core_plugins.text_plugin import TextPlugin


def test_can_be_instantiated():
    assert TextPlugin()


def test_can_be_imported(kernel: Kernel):
    kernel.add_plugin(TextPlugin(), "text_plugin")
    assert not kernel.get_function(plugin_name="text_plugin", function_name="trim").is_prompt


def test_can_be_imported_with_name(kernel: Kernel):
    kernel.add_plugin(TextPlugin(), "text")
    assert not kernel.get_function(plugin_name="text", function_name="trim").is_prompt


def test_can_trim():
    text_plugin = TextPlugin()
    result = text_plugin.trim("  hello world  ")
    assert result == "hello world"


def test_can_trim_start():
    text_plugin = TextPlugin()
    result = text_plugin.trim_start("  hello world  ")
    assert result == "hello world  "


def test_can_trim_end():
    text_plugin = TextPlugin()
    result = text_plugin.trim_end("  hello world  ")
    assert result == "  hello world"


def test_can_lower():
    text_plugin = TextPlugin()
    result = text_plugin.lowercase("  HELLO WORLD  ")
    assert result == "  hello world  "


def test_can_upper():
    text_plugin = TextPlugin()
    result = text_plugin.uppercase("  hello world  ")
    assert result == "  HELLO WORLD  "

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
