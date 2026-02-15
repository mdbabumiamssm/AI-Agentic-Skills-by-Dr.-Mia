# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from unittest.mock import MagicMock

from pytest_mock import MockerFixture

from langchain_perplexity import PerplexitySearchResults


def test_search_tool_run(mocker: MockerFixture) -> None:
    tool = PerplexitySearchResults(pplx_api_key="test")

    mock_result = MagicMock()
    mock_result.title = "Test Title"
    mock_result.url = "http://test.com"
    mock_result.snippet = "Test snippet"
    mock_result.date = "2023-01-01"
    mock_result.last_updated = "2023-01-02"

    mock_response = MagicMock()
    mock_response.results = [mock_result]

    mock_create = MagicMock(return_value=mock_response)
    mocker.patch.object(tool.client.search, "create", mock_create)

    result = tool.invoke("query")

    # result should be a list of dicts (converted by tool) or str if string output
    # By default, tool.invoke returns the output of _run.
    assert isinstance(result, list)
    assert len(result) == 1
    assert result[0]["title"] == "Test Title"

    mock_create.assert_called_once_with(
        query="query",
        max_results=10,
    )

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
