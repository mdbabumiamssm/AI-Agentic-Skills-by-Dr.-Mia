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
    from langchain_community.vectorstores import Redis
    from langchain_community.vectorstores.redis.base import RedisVectorStoreRetriever
    from langchain_community.vectorstores.redis.filters import (
        RedisFilter,
        RedisNum,
        RedisTag,
        RedisText,
    )

# Create a way to dynamically look up deprecated imports.
# Used to consolidate logic for raising deprecation warnings and
# handling optional imports.
DEPRECATED_LOOKUP = {
    "Redis": "langchain_community.vectorstores",
    "RedisFilter": "langchain_community.vectorstores.redis.filters",
    "RedisTag": "langchain_community.vectorstores.redis.filters",
    "RedisText": "langchain_community.vectorstores.redis.filters",
    "RedisNum": "langchain_community.vectorstores.redis.filters",
    "RedisVectorStoreRetriever": "langchain_community.vectorstores.redis.base",
}

_import_attribute = create_importer(__package__, deprecated_lookups=DEPRECATED_LOOKUP)


def __getattr__(name: str) -> Any:
    """Look up attributes dynamically."""
    return _import_attribute(name)


__all__ = [
    "Redis",
    "RedisFilter",
    "RedisNum",
    "RedisTag",
    "RedisText",
    "RedisVectorStoreRetriever",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
