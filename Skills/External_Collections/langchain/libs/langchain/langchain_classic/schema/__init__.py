# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""**Schemas** are the LangChain Base Classes and Interfaces."""

from langchain_core.agents import AgentAction, AgentFinish
from langchain_core.caches import BaseCache
from langchain_core.chat_history import BaseChatMessageHistory
from langchain_core.documents import BaseDocumentTransformer, Document
from langchain_core.exceptions import LangChainException, OutputParserException
from langchain_core.messages import (
    AIMessage,
    BaseMessage,
    ChatMessage,
    FunctionMessage,
    HumanMessage,
    SystemMessage,
    _message_from_dict,
    get_buffer_string,
    messages_from_dict,
    messages_to_dict,
)
from langchain_core.messages.base import message_to_dict
from langchain_core.output_parsers import (
    BaseLLMOutputParser,
    BaseOutputParser,
    StrOutputParser,
)
from langchain_core.outputs import (
    ChatGeneration,
    ChatResult,
    Generation,
    LLMResult,
    RunInfo,
)
from langchain_core.prompt_values import PromptValue
from langchain_core.prompts import BasePromptTemplate, format_document
from langchain_core.retrievers import BaseRetriever
from langchain_core.stores import BaseStore

from langchain_classic.base_memory import BaseMemory

RUN_KEY = "__run"

# Backwards compatibility.
Memory = BaseMemory
_message_to_dict = message_to_dict

__all__ = [
    "RUN_KEY",
    "AIMessage",
    "AgentAction",
    "AgentFinish",
    "BaseCache",
    "BaseChatMessageHistory",
    "BaseDocumentTransformer",
    "BaseLLMOutputParser",
    "BaseMemory",
    "BaseMessage",
    "BaseOutputParser",
    "BasePromptTemplate",
    "BaseRetriever",
    "BaseStore",
    "ChatGeneration",
    "ChatMessage",
    "ChatResult",
    "Document",
    "FunctionMessage",
    "Generation",
    "HumanMessage",
    "LLMResult",
    "LangChainException",
    "Memory",
    "OutputParserException",
    "PromptValue",
    "RunInfo",
    "StrOutputParser",
    "SystemMessage",
    "_message_from_dict",
    "_message_to_dict",
    "format_document",
    "get_buffer_string",
    "message_to_dict",
    "messages_from_dict",
    "messages_to_dict",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
