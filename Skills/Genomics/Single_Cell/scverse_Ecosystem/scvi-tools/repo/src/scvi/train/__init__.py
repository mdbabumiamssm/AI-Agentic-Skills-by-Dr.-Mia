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

from ._callbacks import (
    LoudEarlyStopping,
    SaveCheckpoint,
    ScibCallback,
)
from ._constants import METRIC_KEYS
from ._trainer import Trainer
from ._trainingplans import (
    AdversarialTrainingPlan,
    ClassifierTrainingPlan,
    LowLevelPyroTrainingPlan,
    PyroTrainingPlan,
    SemiSupervisedAdversarialTrainingPlan,
    SemiSupervisedTrainingPlan,
    TrainingPlan,
)
from ._trainrunner import TrainRunner

__all__ = [
    "TrainingPlan",
    "Trainer",
    "PyroTrainingPlan",
    "LowLevelPyroTrainingPlan",
    "SemiSupervisedTrainingPlan",
    "SemiSupervisedAdversarialTrainingPlan",
    "AdversarialTrainingPlan",
    "ClassifierTrainingPlan",
    "TrainRunner",
    "LoudEarlyStopping",
    "SaveCheckpoint",
    "ScibCallback",
    "METRIC_KEYS",
]


def __getattr__(name: str):
    """Lazily provide object. If optional deps are missing, raise a helpful ImportError

    only when object is actually requested.
    """
    if name == "JaxModuleInit":
        error_on_missing_dependencies("flax", "jax", "jaxlib", "optax", "numpyro")
        from ._callbacks import JaxModuleInit as _JaxModuleInit

        return _JaxModuleInit
    if name == "JaxTrainingPlan":
        error_on_missing_dependencies("flax", "jax", "jaxlib", "optax", "numpyro")
        from ._trainingplans import JaxTrainingPlan as _JaxTrainingPlan

        return _JaxTrainingPlan
    raise AttributeError(f"module {__name__!r} has no attribute {name}")

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
