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

from typing import TYPE_CHECKING, Annotated

from semantic_kernel.functions.kernel_function_decorator import kernel_function

if TYPE_CHECKING:
    from semantic_kernel.connectors.search_engine.connector import ConnectorBase


class WebSearchEnginePlugin:
    """A plugin that provides web search engine functionality.

    Usage:
        connector = BingConnector(bing_search_api_key)
        kernel.add_plugin(WebSearchEnginePlugin(connector), plugin_name="WebSearch")

    Examples:
        {{WebSearch.search "What is semantic kernel?"}}
        =>  Returns the first `num_results` number of results for the given search query
            and ignores the first `offset` number of results.
    """

    _connector: "ConnectorBase"

    def __init__(self, connector: "ConnectorBase") -> None:
        """Initializes a new instance of the WebSearchEnginePlugin class."""
        self._connector = connector

    @kernel_function(name="search", description="Performs a web search for a given query")
    async def search(
        self,
        query: Annotated[str, "The search query"],
        num_results: Annotated[int, "The number of search results to return"] = 1,
        offset: Annotated[int, "The number of search results to skip"] = 0,
    ) -> list[str]:
        """Returns the search results of the query provided."""
        return await self._connector.search(query, num_results, offset)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
