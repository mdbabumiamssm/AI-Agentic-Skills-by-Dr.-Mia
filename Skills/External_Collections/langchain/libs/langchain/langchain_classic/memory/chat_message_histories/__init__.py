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
    from langchain_community.chat_message_histories import (
        AstraDBChatMessageHistory,
        CassandraChatMessageHistory,
        ChatMessageHistory,
        CosmosDBChatMessageHistory,
        DynamoDBChatMessageHistory,
        ElasticsearchChatMessageHistory,
        FileChatMessageHistory,
        FirestoreChatMessageHistory,
        MomentoChatMessageHistory,
        MongoDBChatMessageHistory,
        Neo4jChatMessageHistory,
        PostgresChatMessageHistory,
        RedisChatMessageHistory,
        RocksetChatMessageHistory,
        SingleStoreDBChatMessageHistory,
        SQLChatMessageHistory,
        StreamlitChatMessageHistory,
        UpstashRedisChatMessageHistory,
        XataChatMessageHistory,
        ZepChatMessageHistory,
    )

# Create a way to dynamically look up deprecated imports.
# Used to consolidate logic for raising deprecation warnings and
# handling optional imports.
DEPRECATED_LOOKUP = {
    "AstraDBChatMessageHistory": "langchain_community.chat_message_histories",
    "CassandraChatMessageHistory": "langchain_community.chat_message_histories",
    "ChatMessageHistory": "langchain_community.chat_message_histories",
    "CosmosDBChatMessageHistory": "langchain_community.chat_message_histories",
    "DynamoDBChatMessageHistory": "langchain_community.chat_message_histories",
    "ElasticsearchChatMessageHistory": "langchain_community.chat_message_histories",
    "FileChatMessageHistory": "langchain_community.chat_message_histories",
    "FirestoreChatMessageHistory": "langchain_community.chat_message_histories",
    "MomentoChatMessageHistory": "langchain_community.chat_message_histories",
    "MongoDBChatMessageHistory": "langchain_community.chat_message_histories",
    "Neo4jChatMessageHistory": "langchain_community.chat_message_histories",
    "PostgresChatMessageHistory": "langchain_community.chat_message_histories",
    "RedisChatMessageHistory": "langchain_community.chat_message_histories",
    "RocksetChatMessageHistory": "langchain_community.chat_message_histories",
    "SQLChatMessageHistory": "langchain_community.chat_message_histories",
    "SingleStoreDBChatMessageHistory": "langchain_community.chat_message_histories",
    "StreamlitChatMessageHistory": "langchain_community.chat_message_histories",
    "UpstashRedisChatMessageHistory": "langchain_community.chat_message_histories",
    "XataChatMessageHistory": "langchain_community.chat_message_histories",
    "ZepChatMessageHistory": "langchain_community.chat_message_histories",
}

_import_attribute = create_importer(__package__, deprecated_lookups=DEPRECATED_LOOKUP)


def __getattr__(name: str) -> Any:
    """Look up attributes dynamically."""
    return _import_attribute(name)


__all__ = [
    "AstraDBChatMessageHistory",
    "CassandraChatMessageHistory",
    "ChatMessageHistory",
    "CosmosDBChatMessageHistory",
    "DynamoDBChatMessageHistory",
    "ElasticsearchChatMessageHistory",
    "FileChatMessageHistory",
    "FirestoreChatMessageHistory",
    "MomentoChatMessageHistory",
    "MongoDBChatMessageHistory",
    "Neo4jChatMessageHistory",
    "PostgresChatMessageHistory",
    "RedisChatMessageHistory",
    "RocksetChatMessageHistory",
    "SQLChatMessageHistory",
    "SingleStoreDBChatMessageHistory",
    "StreamlitChatMessageHistory",
    "UpstashRedisChatMessageHistory",
    "XataChatMessageHistory",
    "ZepChatMessageHistory",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
