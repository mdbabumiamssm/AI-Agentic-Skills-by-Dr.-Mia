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

import logging

from pytest import mark, raises

from semantic_kernel.exceptions import NamedArgBlockSyntaxError
from semantic_kernel.functions.kernel_arguments import KernelArguments
from semantic_kernel.kernel import Kernel
from semantic_kernel.template_engine.blocks.block_types import BlockTypes
from semantic_kernel.template_engine.blocks.named_arg_block import NamedArgBlock
from semantic_kernel.template_engine.blocks.val_block import ValBlock
from semantic_kernel.template_engine.blocks.var_block import VarBlock

logger = logging.getLogger(__name__)


def test_init_with_var():
    named_arg_block = NamedArgBlock(content="test=$test_var")
    assert named_arg_block.content == "test=$test_var"
    assert named_arg_block.name == "test"
    assert named_arg_block.variable.name == "test_var"
    assert isinstance(named_arg_block.variable, VarBlock)


def test_init_with_val():
    named_arg_block = NamedArgBlock(content="test='test_val'")
    assert named_arg_block.content == "test='test_val'"
    assert named_arg_block.name == "test"
    assert named_arg_block.value.value == "test_val"
    assert isinstance(named_arg_block.value, ValBlock)


def test_type_property():
    named_arg_block = NamedArgBlock(content="test=$test_var")
    assert named_arg_block.type == BlockTypes.NAMED_ARG


@mark.parametrize(
    "content",
    [
        "=$test_var",
        "test=$test-var",
        "test='test_val\"",
        "test=''",
        "test=$",
    ],
    ids=["no_name", "invalid_var", "invalid_val", "empty_val", "empty_var"],
)
def test_syntax_error(content):
    match = content.replace("$", "\\$") if "$" in content else content
    with raises(NamedArgBlockSyntaxError, match=rf".*{match}.*"):
        NamedArgBlock(content=content)


def test_render():
    named_arg_block = NamedArgBlock(content="test=$test_var")
    rendered_value = named_arg_block.render(Kernel(), KernelArguments(test_var="test_value"))
    assert rendered_value == "test_value"


def test_render_variable_not_found():
    named_arg_block = NamedArgBlock(content="test=$test_var")
    rendered_value = named_arg_block.render(Kernel(), KernelArguments())
    assert rendered_value == ""


def test_init_minimal_var():
    block = NamedArgBlock(content="a=$a")
    assert block.name == "a"
    assert block.variable.name == "a"


def test_init_minimal_val():
    block = NamedArgBlock(content="a='a'")
    assert block.name == "a"
    assert block.value.value == "a"


def test_init_empty():
    with raises(NamedArgBlockSyntaxError, match=r".*"):
        NamedArgBlock(content="")


def test_it_trims_spaces():
    assert NamedArgBlock(content="  a=$x  ").content == "a=$x"


def test_it_ignores_spaces_around():
    target = NamedArgBlock(content="  a=$var \n ")
    assert target.content == "a=$var"


def test_it_renders_to_empty_string_without_variables():
    target = NamedArgBlock(content="  a=$var \n ")
    result = target.render(Kernel(), None)
    assert result == ""


def test_it_renders_to_empty_string_if_variable_is_missing():
    target = NamedArgBlock(content="  a=$var \n ")
    result = target.render(Kernel(), KernelArguments(foo="bar"))
    assert result == ""


def test_it_renders_to_variable_value_when_available():
    target = NamedArgBlock(content="  a=$var \n ")
    result = target.render(Kernel(), KernelArguments(foo="bar", var="able"))
    assert result == "able"


def test_it_renders_to_value():
    target = NamedArgBlock(content="  a='var' \n ")
    result = target.render(Kernel(), None)
    assert result == "var"

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
