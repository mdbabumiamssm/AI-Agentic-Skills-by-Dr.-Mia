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

from scvi.data import synthetic_iid
from scvi.external.velovi import VELOVI


@pytest.mark.optional
def test_velovi():
    n_latent = 5
    adata = synthetic_iid()
    adata.layers["spliced"] = adata.X.copy()
    adata.layers["unspliced"] = adata.X.copy()
    VELOVI.setup_anndata(adata, unspliced_layer="unspliced", spliced_layer="spliced")
    model = VELOVI(adata, n_latent=n_latent)
    model.train(1, check_val_every_n_epoch=1, train_size=0.5)
    model.get_latent_representation()
    model.get_velocity()
    model.get_latent_time()
    model.get_state_assignment()
    model.get_expression_fit()
    model.get_directional_uncertainty()
    model.get_permutation_scores(labels_key="labels")

    _ = model.history

    # tests __repr__
    print(model)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
