# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from collections.abc import Generator

from langchain_huggingface.llms import HuggingFacePipeline


def test_huggingface_pipeline_streaming() -> None:
    """Test streaming tokens from huggingface_pipeline."""
    llm = HuggingFacePipeline.from_model_id(
        model_id="openai-community/gpt2",
        task="text-generation",
        pipeline_kwargs={"max_new_tokens": 10},
    )
    generator = llm.stream("Q: How do you say 'hello' in German? A:'", stop=["."])
    stream_results_string = ""
    assert isinstance(generator, Generator)

    for chunk in generator:
        assert isinstance(chunk, str)
        stream_results_string = chunk
    assert len(stream_results_string.strip()) > 0

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
