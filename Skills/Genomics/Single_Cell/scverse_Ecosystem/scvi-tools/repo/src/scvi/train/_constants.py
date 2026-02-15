# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from typing import NamedTuple


class _METRIC_KEYS_NT(NamedTuple):
    TRAINING_KEY: str = "training"
    VALIDATION_KEY: str = "validation"
    # classification
    ACCURACY_KEY: str = "accuracy"
    F1_SCORE_KEY: str = "f1_score"
    CALIBRATION_ERROR_KEY: str = "calibration_error"
    AUROC_KEY: str = "auroc"
    CLASSIFICATION_LOSS_KEY: str = "classification_loss"
    TRUE_LABELS_KEY: str = "true_labels"
    LOGITS_KEY: str = "logits"


METRIC_KEYS = _METRIC_KEYS_NT()

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
