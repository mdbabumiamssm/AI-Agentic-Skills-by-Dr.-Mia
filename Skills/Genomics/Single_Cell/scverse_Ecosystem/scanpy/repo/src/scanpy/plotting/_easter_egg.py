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

from importlib.resources import files

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.image import imread


def dogplot(*_, **__) -> None:
    """Show who’s a good boy."""
    rng = np.random.default_rng()
    n = int(rng.integers(1, 4))
    img_path = files("scanpy.plotting") / f"dogplot_images/doggo_{n}.webp"
    with img_path.open("rb") as f:
        img = imread(f)

    _, ax = plt.subplots(figsize=(3, 3))
    ax.imshow(img)
    ax.set_axis_off()

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
