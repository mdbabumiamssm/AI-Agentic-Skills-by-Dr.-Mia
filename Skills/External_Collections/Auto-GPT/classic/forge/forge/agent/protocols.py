# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from abc import abstractmethod
from typing import TYPE_CHECKING, Awaitable, Generic, Iterator

from forge.models.action import ActionResult, AnyProposal

from .components import AgentComponent

if TYPE_CHECKING:
    from forge.command.command import Command
    from forge.llm.providers import ChatMessage


class DirectiveProvider(AgentComponent):
    def get_constraints(self) -> Iterator[str]:
        return iter([])

    def get_resources(self) -> Iterator[str]:
        return iter([])

    def get_best_practices(self) -> Iterator[str]:
        return iter([])


class CommandProvider(AgentComponent):
    @abstractmethod
    def get_commands(self) -> Iterator["Command"]:
        ...


class MessageProvider(AgentComponent):
    @abstractmethod
    def get_messages(self) -> Iterator["ChatMessage"]:
        ...


class AfterParse(AgentComponent, Generic[AnyProposal]):
    @abstractmethod
    def after_parse(self, result: AnyProposal) -> None | Awaitable[None]:
        ...


class ExecutionFailure(AgentComponent):
    @abstractmethod
    def execution_failure(self, error: Exception) -> None | Awaitable[None]:
        ...


class AfterExecute(AgentComponent):
    @abstractmethod
    def after_execute(self, result: "ActionResult") -> None | Awaitable[None]:
        ...

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
