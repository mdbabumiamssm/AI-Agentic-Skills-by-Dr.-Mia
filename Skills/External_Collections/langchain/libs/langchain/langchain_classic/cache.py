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
    from langchain_community.cache import (
        AstraDBCache,
        AstraDBSemanticCache,
        AzureCosmosDBSemanticCache,
        CassandraCache,
        CassandraSemanticCache,
        FullLLMCache,
        FullMd5LLMCache,
        GPTCache,
        InMemoryCache,
        MomentoCache,
        RedisCache,
        RedisSemanticCache,
        SQLAlchemyCache,
        SQLAlchemyMd5Cache,
        SQLiteCache,
        UpstashRedisCache,
    )

# Create a way to dynamically look up deprecated imports.
# Used to consolidate logic for raising deprecation warnings and
# handling optional imports.
DEPRECATED_LOOKUP = {
    "FullLLMCache": "langchain_community.cache",
    "SQLAlchemyCache": "langchain_community.cache",
    "SQLiteCache": "langchain_community.cache",
    "UpstashRedisCache": "langchain_community.cache",
    "RedisCache": "langchain_community.cache",
    "RedisSemanticCache": "langchain_community.cache",
    "GPTCache": "langchain_community.cache",
    "MomentoCache": "langchain_community.cache",
    "InMemoryCache": "langchain_community.cache",
    "CassandraCache": "langchain_community.cache",
    "CassandraSemanticCache": "langchain_community.cache",
    "FullMd5LLMCache": "langchain_community.cache",
    "SQLAlchemyMd5Cache": "langchain_community.cache",
    "AstraDBCache": "langchain_community.cache",
    "AstraDBSemanticCache": "langchain_community.cache",
    "AzureCosmosDBSemanticCache": "langchain_community.cache",
}

_import_attribute = create_importer(__package__, deprecated_lookups=DEPRECATED_LOOKUP)


def __getattr__(name: str) -> Any:
    """Look up attributes dynamically."""
    return _import_attribute(name)


__all__ = [
    "AstraDBCache",
    "AstraDBSemanticCache",
    "AzureCosmosDBSemanticCache",
    "CassandraCache",
    "CassandraSemanticCache",
    "FullLLMCache",
    "FullMd5LLMCache",
    "GPTCache",
    "InMemoryCache",
    "MomentoCache",
    "RedisCache",
    "RedisSemanticCache",
    "SQLAlchemyCache",
    "SQLAlchemyMd5Cache",
    "SQLiteCache",
    "UpstashRedisCache",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
