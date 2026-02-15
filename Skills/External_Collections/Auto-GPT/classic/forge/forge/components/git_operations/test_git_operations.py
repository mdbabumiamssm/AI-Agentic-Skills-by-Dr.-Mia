# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import pytest
from git.exc import GitCommandError
from git.repo.base import Repo

from forge.file_storage.base import FileStorage
from forge.utils.exceptions import CommandExecutionError

from . import GitOperationsComponent


@pytest.fixture
def mock_clone_from(mocker):
    return mocker.patch.object(Repo, "clone_from")


@pytest.fixture
def git_ops_component():
    return GitOperationsComponent()


def test_clone_auto_gpt_repository(
    git_ops_component: GitOperationsComponent,
    storage: FileStorage,
    mock_clone_from,
):
    mock_clone_from.return_value = None

    repo = "github.com/Significant-Gravitas/Auto-GPT.git"
    scheme = "https://"
    url = scheme + repo
    clone_path = storage.get_path("auto-gpt-repo")

    expected_output = f"Cloned {url} to {clone_path}"

    clone_result = git_ops_component.clone_repository(url, clone_path)

    assert clone_result == expected_output
    mock_clone_from.assert_called_once_with(
        url=f"{scheme}{git_ops_component.config.github_username}:{git_ops_component.config.github_api_key}@{repo}",  # noqa: E501
        to_path=clone_path,
    )


def test_clone_repository_error(
    git_ops_component: GitOperationsComponent,
    storage: FileStorage,
    mock_clone_from,
):
    url = "https://github.com/this-repository/does-not-exist.git"
    clone_path = storage.get_path("does-not-exist")

    mock_clone_from.side_effect = GitCommandError(
        "clone", "fatal: repository not found", ""
    )

    with pytest.raises(CommandExecutionError):
        git_ops_component.clone_repository(url, clone_path)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
