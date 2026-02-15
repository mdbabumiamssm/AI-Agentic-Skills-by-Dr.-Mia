# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from langchain_core.prompts import PromptTemplate

template = """You are a teacher coming up with questions to ask on a quiz.
Given the following document, please generate a question and answer based on that document.

Example Format:
<Begin Document>
...
<End Document>
QUESTION: question here
ANSWER: answer here

These questions should be detailed and be based explicitly on information in the document. Begin!

<Begin Document>
{doc}
<End Document>"""  # noqa: E501
PROMPT = PromptTemplate(
    input_variables=["doc"],
    template=template,
)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
