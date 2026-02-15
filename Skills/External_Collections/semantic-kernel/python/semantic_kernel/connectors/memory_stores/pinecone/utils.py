# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

# Copyright (c) Microsoft. All rights reserved.

import numpy
from pinecone import Vector

from semantic_kernel.memory.memory_record import MemoryRecord


def build_payload(record: MemoryRecord) -> dict:
    """Builds a metadata payload to be sent to Pinecone from a MemoryRecord."""
    payload: dict = {}
    if record._text:
        payload["text"] = record._text
    if record._description:
        payload["description"] = record._description
    if record._additional_metadata:
        payload["additional_metadata"] = record._additional_metadata
    return payload


def parse_payload(record: Vector, with_embeddings: bool) -> MemoryRecord:
    """Parses a record from Pinecone into a MemoryRecord."""
    payload = record.metadata
    description = payload.get("description", None)
    text = payload.get("text", None)
    additional_metadata = payload.get("additional_metadata", None)
    return MemoryRecord.local_record(
        id=record.id,
        description=description,
        text=text,
        additional_metadata=additional_metadata,
        embedding=record.values if with_embeddings else numpy.array([]),
    )

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
