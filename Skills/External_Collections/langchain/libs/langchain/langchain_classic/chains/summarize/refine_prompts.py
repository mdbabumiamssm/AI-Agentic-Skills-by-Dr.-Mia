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

REFINE_PROMPT_TMPL = """\
Your job is to produce a final summary.
We have provided an existing summary up to a certain point: {existing_answer}
We have the opportunity to refine the existing summary (only if needed) with some more context below.
------------
{text}
------------
Given the new context, refine the original summary.
If the context isn't useful, return the original summary.\
"""  # noqa: E501
REFINE_PROMPT = PromptTemplate.from_template(REFINE_PROMPT_TMPL)


prompt_template = """Write a concise summary of the following:


"{text}"


CONCISE SUMMARY:"""
PROMPT = PromptTemplate.from_template(prompt_template)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
