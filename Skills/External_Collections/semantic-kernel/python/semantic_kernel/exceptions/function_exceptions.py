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
from semantic_kernel.exceptions.kernel_exceptions import KernelException


class FunctionException(KernelException):
    """Base class for all function exceptions."""

    pass


class FunctionSyntaxError(FunctionException):
    """Base class for all function syntax exceptions."""

    pass


class FunctionInitializationError(FunctionException):
    """An error occurred while initializing the function."""

    def __init__(self, message: str):
        """Adds the context of the error to the generic message."""
        super().__init__("KernelFunction failed to initialize: " + message)


class PluginInitializationError(FunctionException):
    """An error occurred while initializing the plugin."""

    pass


class PluginInvalidNameError(FunctionSyntaxError):
    """An error occurred while validating the plugin name."""

    pass


class FunctionInvalidNameError(FunctionSyntaxError):
    """An error occurred while validating the function name."""

    pass


class FunctionInvalidParamNameError(FunctionSyntaxError):
    """An error occurred while validating the function parameter name."""

    pass


class FunctionNameNotUniqueError(FunctionSyntaxError):
    """An error occurred while validating the function name."""

    pass


class FunctionExecutionException(FunctionException):
    """Base class for all function execution exceptions."""

    pass


class FunctionResultError(FunctionException):
    """An error occurred while validating the function result."""

    pass


class FunctionInvalidParameterConfiguration(FunctionException):
    """An error occurred while validating the function parameter configuration."""

    pass


class PromptRenderingException(FunctionException):
    """An error occurred while rendering a prompt."""

    pass


__all__ = [
    "FunctionException",
    "FunctionExecutionException",
    "FunctionInitializationError",
    "FunctionInvalidNameError",
    "FunctionInvalidParamNameError",
    "FunctionNameNotUniqueError",
    "FunctionResultError",
    "FunctionSyntaxError",
    "PluginInitializationError",
    "PluginInvalidNameError",
    "PromptRenderingException",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
