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


class _SCVI_HUB_NT(NamedTuple):
    HF_LIBRARY_NAME: str = "scvi-tools"
    MAX_HF_UPLOAD_SIZE: int = 5e9  # 5GB

    # file names
    METADATA_FILE_NAME: str = "_scvi_required_metadata.json"
    MODEL_CARD_FILE_NAME: str = "README.md"

    # model card defaults
    DEFAULT_MISSING_FIELD: str = "To be added..."
    DEFAULT_NA_FIELD: str = "Not provided by uploader"
    DEFAULT_PARENT_MODULE: str = "scvi.model"

    # model card tags
    MODEL_CLS_NAME_TAG: str = "model_cls_name:{}"
    SCVI_VERSION_TAG: str = "scvi_version:{}"
    ANNDATA_VERSION_TAG: str = "anndata_version:{}"
    MODALITY_TAG: str = "modality:{}"
    TISSUE_TAG: str = "tissue:{}"
    ANNOTATED_TAG: str = "annotated:{}"


_SCVI_HUB = _SCVI_HUB_NT()

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
