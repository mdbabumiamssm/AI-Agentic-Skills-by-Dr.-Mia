# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from langchain_core.runnables import __all__

EXPECTED_ALL = [
    "chain",
    "AddableDict",
    "ConfigurableField",
    "ConfigurableFieldSingleOption",
    "ConfigurableFieldMultiOption",
    "ConfigurableFieldSpec",
    "ensure_config",
    "run_in_executor",
    "patch_config",
    "RouterInput",
    "RouterRunnable",
    "Runnable",
    "RunnableSerializable",
    "RunnableBinding",
    "RunnableBranch",
    "RunnableConfig",
    "RunnableGenerator",
    "RunnableLambda",
    "RunnableMap",
    "RunnableParallel",
    "RunnablePassthrough",
    "RunnableAssign",
    "RunnablePick",
    "RunnableSequence",
    "RunnableWithFallbacks",
    "RunnableWithMessageHistory",
    "get_config_list",
    "aadd",
    "add",
]


def test_all_imports() -> None:
    assert set(__all__) == set(EXPECTED_ALL)


def test_imports_for_specific_funcs() -> None:
    """Test that a few specific imports in more internal namespaces."""
    # create_model implementation has been moved to langchain_core.utils.pydantic
    from langchain_core.runnables.utils import (  # type: ignore[attr-defined] # noqa: F401,PLC0415
        create_model,
    )

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
