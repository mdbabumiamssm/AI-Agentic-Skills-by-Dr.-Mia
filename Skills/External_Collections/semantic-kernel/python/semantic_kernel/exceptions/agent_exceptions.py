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


class AgentException(KernelException):
    """Base class for all agent exceptions."""

    pass


class AgentFileNotFoundException(AgentException):
    """The requested file was not found."""

    pass


class AgentInitializationException(AgentException):
    """An error occurred while initializing the agent."""

    pass


class AgentExecutionException(AgentException):
    """An error occurred while executing the agent."""

    pass


class AgentInvokeException(AgentException):
    """An error occurred while invoking the agent."""

    pass


class AgentChatException(AgentException):
    """An error occurred while invoking the agent chat."""

    pass


class AgentChatHistoryReducerException(AgentException):
    """An error occurred while reducing the chat history."""

    pass


class AgentThreadInitializationException(AgentException):
    """An error occurred while initializing the agent thread."""

    pass


class AgentThreadOperationException(AgentException):
    """An error occurred while performing an operation on the agent thread."""

    pass

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
