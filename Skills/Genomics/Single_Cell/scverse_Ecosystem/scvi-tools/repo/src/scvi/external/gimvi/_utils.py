# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import os
import pickle
from typing import Literal

import numpy as np
import torch
from anndata import AnnData, read_h5ad

from scvi.data._download import _download


def _load_legacy_saved_gimvi_files(
    dir_path: str,
    file_name_prefix: str,
    load_seq_adata: bool,
    load_spatial_adata: bool,
) -> tuple[dict, np.ndarray, np.ndarray, dict, AnnData | None, AnnData | None]:
    model_path = os.path.join(dir_path, f"{file_name_prefix}model_params.pt")
    setup_dict_path = os.path.join(dir_path, f"{file_name_prefix}attr.pkl")
    seq_var_names_path = os.path.join(dir_path, f"{file_name_prefix}var_names_seq.csv")
    spatial_var_names_path = os.path.join(dir_path, f"{file_name_prefix}var_names_spatial.csv")

    model_state_dict = torch.load(model_path, map_location="cpu", weights_only=False)

    seq_var_names = np.genfromtxt(seq_var_names_path, delimiter=",", dtype=str)
    spatial_var_names = np.genfromtxt(spatial_var_names_path, delimiter=",", dtype=str)

    with open(setup_dict_path, "rb") as handle:
        attr_dict = pickle.load(handle)

    adata_seq, adata_spatial = None, None
    if load_seq_adata:
        seq_data_path = os.path.join(dir_path, f"{file_name_prefix}adata_seq.h5ad")
        if os.path.exists(seq_data_path):
            adata_seq = read_h5ad(seq_data_path)
        elif not os.path.exists(seq_data_path):
            Warning(
                "Save path contains no saved anndata and no adata was passed. "
                "Model will be loaded without anndata."
            )
    if load_spatial_adata:
        spatial_data_path = os.path.join(dir_path, f"{file_name_prefix}adata_spatial.h5ad")
        if os.path.exists(spatial_data_path):
            adata_spatial = read_h5ad(spatial_data_path)
        elif not os.path.exists(spatial_data_path):
            Warning(
                "Save path contains no saved anndata and no adata was passed. "
                "Model will be loaded without anndata."
            )

    return (
        model_state_dict,
        seq_var_names,
        spatial_var_names,
        attr_dict,
        adata_seq,
        adata_spatial,
    )


def _load_saved_gimvi_files(
    dir_path: str,
    load_seq_adata: bool,
    load_spatial_adata: bool,
    prefix: str | None = None,
    map_location: Literal["cpu", "cuda"] | None = None,
    backup_url: str | None = None,
) -> tuple[dict, dict, np.ndarray, np.ndarray, dict, AnnData | None, AnnData | None]:
    file_name_prefix = prefix or ""

    model_file_name = f"{file_name_prefix}model.pt"
    model_path = os.path.join(dir_path, model_file_name)
    try:
        _download(backup_url, dir_path, model_file_name)
        model = torch.load(model_path, map_location=map_location, weights_only=False)
    except FileNotFoundError as exc:
        raise ValueError(
            f"Failed to load model file at {model_path}. "
            "If attempting to load a saved model from <v0.15.0, please use the util function "
            "`convert_legacy_save` to convert to an updated format."
        ) from exc

    model_state_dict = model["model_state_dict"]
    seq_var_names = model["seq_var_names"]
    spatial_var_names = model["spatial_var_names"]
    attr_dict = model["attr_dict"]

    adata_seq, adata_spatial = None, None
    if load_seq_adata:
        seq_data_path = os.path.join(dir_path, f"{file_name_prefix}adata_seq.h5ad")
        if os.path.exists(seq_data_path):
            adata_seq = read_h5ad(seq_data_path)
        elif not os.path.exists(seq_data_path):
            Warning(
                "Save path contains no saved anndata and no adata was passed. "
                "Model will be loaded without anndata."
            )
    if load_spatial_adata:
        spatial_data_path = os.path.join(dir_path, f"{file_name_prefix}adata_spatial.h5ad")
        if os.path.exists(spatial_data_path):
            adata_spatial = read_h5ad(spatial_data_path)
        elif not os.path.exists(spatial_data_path):
            Warning(
                "Save path contains no saved anndata and no adata was passed. "
                "Model will be loaded without anndata."
            )

    return (
        attr_dict,
        seq_var_names,
        spatial_var_names,
        model_state_dict,
        adata_seq,
        adata_spatial,
    )

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
