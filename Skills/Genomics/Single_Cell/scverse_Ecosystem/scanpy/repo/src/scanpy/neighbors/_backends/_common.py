# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from __future__ import annotations


class TransformerChecksMixin:
    def _transform_checks(self, x, /, *fitted_props, **check_params):
        from sklearn.utils.validation import check_is_fitted

        if x is not None:
            x = self._validate_data(x, reset=False, **check_params)
        check_is_fitted(self, *fitted_props)
        return x

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
