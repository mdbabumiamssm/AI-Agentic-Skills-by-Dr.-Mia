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

DEFAULT_REFINE_PROMPT_TMPL = (
    "The original question is as follows: {question}\n"
    "We have provided an existing answer, including sources: {existing_answer}\n"
    "We have the opportunity to refine the existing answer"
    "(only if needed) with some more context below.\n"
    "------------\n"
    "{context_str}\n"
    "------------\n"
    "Given the new context, refine the original answer to better "
    "answer the question. "
    "If you do update it, please update the sources as well. "
    "If the context isn't useful, return the original answer."
)
DEFAULT_REFINE_PROMPT = PromptTemplate(
    input_variables=["question", "existing_answer", "context_str"],
    template=DEFAULT_REFINE_PROMPT_TMPL,
)


DEFAULT_TEXT_QA_PROMPT_TMPL = (
    "Context information is below. \n"
    "---------------------\n"
    "{context_str}"
    "\n---------------------\n"
    "Given the context information and not prior knowledge, "
    "answer the question: {question}\n"
)
DEFAULT_TEXT_QA_PROMPT = PromptTemplate(
    input_variables=["context_str", "question"], template=DEFAULT_TEXT_QA_PROMPT_TMPL
)

EXAMPLE_PROMPT = PromptTemplate(
    template="Content: {page_content}\nSource: {source}",
    input_variables=["page_content", "source"],
)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
