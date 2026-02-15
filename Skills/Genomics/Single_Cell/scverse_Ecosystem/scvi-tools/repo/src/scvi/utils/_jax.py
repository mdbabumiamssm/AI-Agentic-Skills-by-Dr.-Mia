# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from collections.abc import Callable

from scvi.utils import dependencies


@dependencies("jax")
def device_selecting_PRNGKey(use_cpu: bool = True) -> Callable:
    """Returns a PRNGKey that is either on CPU or GPU."""
    # if key is generated on CPU, model params will be on CPU
    import jax
    from jax import random

    if use_cpu is True:

        def key(i: int):
            return jax.device_put(random.PRNGKey(i), jax.devices("cpu")[0])
    else:
        # dummy function
        def key(i: int):
            return random.PRNGKey(i)

    return key

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
