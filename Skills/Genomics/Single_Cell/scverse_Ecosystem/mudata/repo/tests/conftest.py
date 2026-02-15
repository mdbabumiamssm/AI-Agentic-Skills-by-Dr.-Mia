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


@pytest.fixture(scope="module")
def filepath_h5mu(tmpdir_factory):
    return str(tmpdir_factory.mktemp("tmp_test_dir").join("testA.h5mu"))


@pytest.fixture(scope="module")
def filepath2_h5mu(tmpdir_factory):
    return str(tmpdir_factory.mktemp("tmp_test_dir").join("testB.h5mu"))


@pytest.fixture(scope="module")
def filepath_hdf5(tmpdir_factory):
    return str(tmpdir_factory.mktemp("tmp_mofa_dir").join("mofa_pytest.hdf5"))


@pytest.fixture(scope="module")
def filepath_zarr(tmpdir_factory):
    return str(tmpdir_factory.mktemp("tmp_test_dir").join("testA.zarr"))


@pytest.fixture(scope="module")
def filepath2_zarr(tmpdir_factory):
    return str(tmpdir_factory.mktemp("tmp_test_dir").join("testB.zarr"))

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
