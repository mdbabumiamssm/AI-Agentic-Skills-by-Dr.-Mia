# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from scvi.utils import error_on_missing_dependencies

from ._beta_binomial import BetaBinomial
from ._negative_binomial import (
    NegativeBinomial,
    NegativeBinomialMixture,
    Poisson,
    ZeroInflatedNegativeBinomial,
)
from ._normal import Normal

__all__ = [
    "NegativeBinomial",
    "NegativeBinomialMixture",
    "ZeroInflatedNegativeBinomial",
    "Poisson",
    "BetaBinomial",
    "Normal",
]


def __getattr__(name: str):
    """
    Lazily provide object. If optional deps are missing, raise a helpful ImportError

    only when object is actually requested.
    """
    if name == "JaxNegativeBinomialMeanDisp":
        error_on_missing_dependencies("jax", "numpyro")
        from ._negative_binomial import JaxNegativeBinomialMeanDisp as _JaxNegativeBinomialMeanDisp

        return _JaxNegativeBinomialMeanDisp
    raise AttributeError(f"module {__name__!r} has no attribute {name}")

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
