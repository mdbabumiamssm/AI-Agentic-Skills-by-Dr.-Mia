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
# <defineClass>
import math
from typing import Annotated

from semantic_kernel.functions.kernel_function_decorator import kernel_function


class Math:
    # </defineClass>
    """
    Description: MathPlugin provides a set of functions to make Math calculations.

    Usage:
        kernel.add_plugin(MathPlugin(), plugin_name="math")

    Examples:
        {{math.Add}} => Returns the sum of input and amount (provided in the KernelArguments)
        {{math.Subtract}} => Returns the difference of input and amount (provided in the KernelArguments)
        {{math.Multiply}} => Returns the multiplication of input and number2 (provided in the KernelArguments)
        {{math.Divide}} => Returns the division of input and number2 (provided in the KernelArguments)
    """

    @kernel_function(
        description="Divide two numbers.",
        name="Divide",
    )
    def divide(
        self,
        number1: Annotated[float, "the first number to divide from"],
        number2: Annotated[float, "the second number to by"],
    ) -> Annotated[float, "The output is a float"]:
        return float(number1) / float(number2)

    @kernel_function(
        description="Multiply two numbers. When increasing by a percentage, don't forget to add 1 to the percentage.",
        name="Multiply",
    )
    def multiply(
        self,
        number1: Annotated[float, "the first number to multiply"],
        number2: Annotated[float, "the second number to multiply"],
    ) -> Annotated[float, "The output is a float"]:
        return float(number1) * float(number2)

    # <defineFunction>
    @kernel_function(
        description="Takes the square root of a number",
        name="Sqrt",
    )
    def square_root(
        self,
        number1: Annotated[float, "the number to take the square root of"],
    ) -> Annotated[float, "The output is a float"]:
        return math.sqrt(float(number1))

    # </defineFunction>

    @kernel_function(name="Add")
    def add(
        self,
        number1: Annotated[float, "the first number to add"],
        number2: Annotated[float, "the second number to add"],
    ) -> Annotated[float, "the output is a float"]:
        return float(number1) + float(number2)

    @kernel_function(
        description="Subtracts value to a value",
        name="Subtract",
    )
    def subtract(
        self,
        number1: Annotated[float, "the first number"],
        number2: Annotated[float, "the number to subtract"],
    ) -> Annotated[float, "the output is a float"]:
        return float(number1) - float(number2)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
