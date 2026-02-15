# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from typing import TYPE_CHECKING, Any

from langchain_classic._api import create_importer

if TYPE_CHECKING:
    from langchain_community.vectorstores.redis.schema import (
        FlatVectorField,
        HNSWVectorField,
        NumericFieldSchema,
        RedisDistanceMetric,
        RedisField,
        RedisModel,
        RedisVectorField,
        TagFieldSchema,
        TextFieldSchema,
        read_schema,
    )

# Create a way to dynamically look up deprecated imports.
# Used to consolidate logic for raising deprecation warnings and
# handling optional imports.
DEPRECATED_LOOKUP = {
    "RedisDistanceMetric": "langchain_community.vectorstores.redis.schema",
    "RedisField": "langchain_community.vectorstores.redis.schema",
    "TextFieldSchema": "langchain_community.vectorstores.redis.schema",
    "TagFieldSchema": "langchain_community.vectorstores.redis.schema",
    "NumericFieldSchema": "langchain_community.vectorstores.redis.schema",
    "RedisVectorField": "langchain_community.vectorstores.redis.schema",
    "FlatVectorField": "langchain_community.vectorstores.redis.schema",
    "HNSWVectorField": "langchain_community.vectorstores.redis.schema",
    "RedisModel": "langchain_community.vectorstores.redis.schema",
    "read_schema": "langchain_community.vectorstores.redis.schema",
}

_import_attribute = create_importer(__package__, deprecated_lookups=DEPRECATED_LOOKUP)


def __getattr__(name: str) -> Any:
    """Look up attributes dynamically."""
    return _import_attribute(name)


__all__ = [
    "FlatVectorField",
    "HNSWVectorField",
    "NumericFieldSchema",
    "RedisDistanceMetric",
    "RedisField",
    "RedisModel",
    "RedisVectorField",
    "TagFieldSchema",
    "TextFieldSchema",
    "read_schema",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
