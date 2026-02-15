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

from typing import Literal

import anndata
import mudata
import torch

from scvi.utils import is_package_installed

Number = int | float
AnnOrMuData = anndata.AnnData | mudata.MuData
if is_package_installed("jax"):
    import jax.numpy as jnp

    Tensor = torch.Tensor | jnp.ndarray
else:
    Tensor = torch.Tensor
LossRecord = dict[str, Tensor] | Tensor
MinifiedDataType = Literal["latent_posterior_parameters"]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
