# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from langchain_prompty.core import InvokerFactory
from langchain_prompty.langchain import create_chat_prompt
from langchain_prompty.parsers import PromptyChatParser
from langchain_prompty.renderers import MustacheRenderer

InvokerFactory().register_renderer("mustache", MustacheRenderer)
InvokerFactory().register_parser("prompty.chat", PromptyChatParser)

__all__ = ["create_chat_prompt"]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
