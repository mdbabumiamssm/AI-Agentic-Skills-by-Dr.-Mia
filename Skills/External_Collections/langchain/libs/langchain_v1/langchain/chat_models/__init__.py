# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Entrypoint to using [chat models](https://docs.langchain.com/oss/python/langchain/models) in LangChain."""  # noqa: E501

from langchain_core.language_models import BaseChatModel

from langchain.chat_models.base import init_chat_model

__all__ = ["BaseChatModel", "init_chat_model"]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
