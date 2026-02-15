# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from collections.abc import Sequence
from inspect import signature

from langchain_core.callbacks import Callbacks
from langchain_core.documents import (
    BaseDocumentCompressor,
    BaseDocumentTransformer,
    Document,
)
from pydantic import ConfigDict


class DocumentCompressorPipeline(BaseDocumentCompressor):
    """Document compressor that uses a pipeline of Transformers."""

    transformers: list[BaseDocumentTransformer | BaseDocumentCompressor]
    """List of document filters that are chained together and run in sequence."""

    model_config = ConfigDict(
        arbitrary_types_allowed=True,
    )

    def compress_documents(
        self,
        documents: Sequence[Document],
        query: str,
        callbacks: Callbacks | None = None,
    ) -> Sequence[Document]:
        """Transform a list of documents."""
        for _transformer in self.transformers:
            if isinstance(_transformer, BaseDocumentCompressor):
                accepts_callbacks = (
                    signature(_transformer.compress_documents).parameters.get(
                        "callbacks",
                    )
                    is not None
                )
                if accepts_callbacks:
                    documents = _transformer.compress_documents(
                        documents,
                        query,
                        callbacks=callbacks,
                    )
                else:
                    documents = _transformer.compress_documents(documents, query)
            elif isinstance(_transformer, BaseDocumentTransformer):
                documents = _transformer.transform_documents(documents)
            else:
                msg = f"Got unexpected transformer type: {_transformer}"  # type: ignore[unreachable]
                raise ValueError(msg)  # noqa: TRY004
        return documents

    async def acompress_documents(
        self,
        documents: Sequence[Document],
        query: str,
        callbacks: Callbacks | None = None,
    ) -> Sequence[Document]:
        """Compress retrieved documents given the query context."""
        for _transformer in self.transformers:
            if isinstance(_transformer, BaseDocumentCompressor):
                accepts_callbacks = (
                    signature(_transformer.acompress_documents).parameters.get(
                        "callbacks",
                    )
                    is not None
                )
                if accepts_callbacks:
                    documents = await _transformer.acompress_documents(
                        documents,
                        query,
                        callbacks=callbacks,
                    )
                else:
                    documents = await _transformer.acompress_documents(documents, query)
            elif isinstance(_transformer, BaseDocumentTransformer):
                documents = await _transformer.atransform_documents(documents)
            else:
                msg = f"Got unexpected transformer type: {_transformer}"  # type: ignore[unreachable]
                raise ValueError(msg)  # noqa: TRY004
        return documents

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
