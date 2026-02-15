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
import tarfile

import pytest
import scanpy as sc

from scvi.data import dataset_10x, organize_cite_seq_10x

from .utils import unsupervised_training_one_epoch


def test_populate_and_train_one_v1(save_path):
    sp = os.path.join(save_path, "10X")
    dataset = dataset_10x(
        dataset_name="cd4_t_helper",
        remove_extracted_data=True,
        save_path=sp,
    )
    unsupervised_training_one_epoch(dataset)


def test_brain_small(save_path):
    sp = os.path.join(save_path, "10X")
    dataset = dataset_10x(
        dataset_name="neuron_9k",
        save_path=sp,
        remove_extracted_data=True,
    )
    unsupervised_training_one_epoch(dataset)


def test_pbmc_cite(save_path):
    file_path = os.path.join(
        save_path, "10X/pbmc_10k_protein_v3/filtered_feature_bc_matrix.tar.gz"
    )
    sp = os.path.join(save_path, "10X/pbmc_10k_protein_v3/")
    tar = tarfile.open(file_path, "r:gz")
    tar.extractall(path=sp)
    tar.close()
    dataset = sc.read_10x_mtx(os.path.join(sp, "filtered_feature_bc_matrix"), gex_only=False)
    organize_cite_seq_10x(dataset)
    unsupervised_training_one_epoch(dataset)


@pytest.mark.internet
def test_download_dataset_10x(save_path):
    dataset_10x("hgmm_1k_v3", save_path=save_path)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
