# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from typing import Any


def get_prompt_input_key(inputs: dict[str, Any], memory_variables: list[str]) -> str:
    """Get the prompt input key.

    Args:
        inputs: Dict[str, Any]
        memory_variables: List[str]

    Returns:
        A prompt input key.
    """
    # "stop" is a special key that can be passed as input but is not used to
    # format the prompt.
    prompt_input_keys = list(set(inputs).difference([*memory_variables, "stop"]))
    if len(prompt_input_keys) != 1:
        msg = f"One input key expected got {prompt_input_keys}"
        raise ValueError(msg)
    return prompt_input_keys[0]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
