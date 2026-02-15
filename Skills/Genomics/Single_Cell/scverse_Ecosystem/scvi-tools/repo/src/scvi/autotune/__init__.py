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

error_on_missing_dependencies("hyperopt", "ray.tune")

from ._experiment import AutotuneExperiment, ScibTuneReportCheckpointCallback  # noqa: E402
from ._tune import run_autotune  # noqa: E402

__all__ = ["AutotuneExperiment", "ScibTuneReportCheckpointCallback", "run_autotune"]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
