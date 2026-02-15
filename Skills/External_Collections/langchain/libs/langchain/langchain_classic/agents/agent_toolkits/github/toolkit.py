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
    from langchain_community.agent_toolkits.github.toolkit import (
        BranchName,
        CommentOnIssue,
        CreateFile,
        CreatePR,
        CreateReviewRequest,
        DeleteFile,
        DirectoryPath,
        GetIssue,
        GetPR,
        GitHubToolkit,
        NoInput,
        ReadFile,
        SearchCode,
        SearchIssuesAndPRs,
        UpdateFile,
    )

# Create a way to dynamically look up deprecated imports.
# Used to consolidate logic for raising deprecation warnings and
# handling optional imports.
DEPRECATED_LOOKUP = {
    "NoInput": "langchain_community.agent_toolkits.github.toolkit",
    "GetIssue": "langchain_community.agent_toolkits.github.toolkit",
    "CommentOnIssue": "langchain_community.agent_toolkits.github.toolkit",
    "GetPR": "langchain_community.agent_toolkits.github.toolkit",
    "CreatePR": "langchain_community.agent_toolkits.github.toolkit",
    "CreateFile": "langchain_community.agent_toolkits.github.toolkit",
    "ReadFile": "langchain_community.agent_toolkits.github.toolkit",
    "UpdateFile": "langchain_community.agent_toolkits.github.toolkit",
    "DeleteFile": "langchain_community.agent_toolkits.github.toolkit",
    "DirectoryPath": "langchain_community.agent_toolkits.github.toolkit",
    "BranchName": "langchain_community.agent_toolkits.github.toolkit",
    "SearchCode": "langchain_community.agent_toolkits.github.toolkit",
    "CreateReviewRequest": "langchain_community.agent_toolkits.github.toolkit",
    "SearchIssuesAndPRs": "langchain_community.agent_toolkits.github.toolkit",
    "GitHubToolkit": "langchain_community.agent_toolkits.github.toolkit",
}

_import_attribute = create_importer(__package__, deprecated_lookups=DEPRECATED_LOOKUP)


def __getattr__(name: str) -> Any:
    """Look up attributes dynamically."""
    return _import_attribute(name)


__all__ = [
    "BranchName",
    "CommentOnIssue",
    "CreateFile",
    "CreatePR",
    "CreateReviewRequest",
    "DeleteFile",
    "DirectoryPath",
    "GetIssue",
    "GetPR",
    "GitHubToolkit",
    "NoInput",
    "ReadFile",
    "SearchCode",
    "SearchIssuesAndPRs",
    "UpdateFile",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
