# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

PREFIX = """You are an agent designed to answer questions about sets of documents.
You have access to tools for interacting with the documents, and the inputs to the tools are questions.
Sometimes, you will be asked to provide sources for your questions, in which case you should use the appropriate tool to do so.
If the question does not seem relevant to any of the tools provided, just return "I don't know" as the answer.
"""  # noqa: E501

ROUTER_PREFIX = """You are an agent designed to answer questions.
You have access to tools for interacting with different sources, and the inputs to the tools are questions.
Your main task is to decide which of the tools is relevant for answering question at hand.
For complex questions, you can break the question down into sub questions and use tools to answers the sub questions.
"""  # noqa: E501

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
