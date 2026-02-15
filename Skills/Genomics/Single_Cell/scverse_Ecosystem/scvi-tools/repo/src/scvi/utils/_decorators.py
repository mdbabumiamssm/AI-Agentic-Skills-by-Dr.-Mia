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
from functools import wraps

from scvi.data._constants import ADATA_MINIFY_TYPE


def unsupported_if_adata_minified(fn: Callable) -> Callable:
    """Decorator to raise an error if the model's `adata` is minified."""

    @wraps(fn)
    def wrapper(self, *args, **kwargs):
        if getattr(self, "minified_data_type", None) == ADATA_MINIFY_TYPE.LATENT_POSTERIOR:
            raise ValueError(
                f"The {fn.__qualname__} function currently does not support minified data."
            )
        return fn(self, *args, **kwargs)

    return wrapper

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
