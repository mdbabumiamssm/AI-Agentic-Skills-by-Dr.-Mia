# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

# for backwards compatibility, this was moved to scvi.data
from scvi.data import AnnTorchDataset

from ._ann_dataloader import AnnDataLoader
from ._anncollection import CollectionAdapter
from ._concat_dataloader import ConcatDataLoader
from ._custom_dataloaders import MappedCollectionDataModule, TileDBDataModule
from ._data_splitting import (
    DataSplitter,
    DeviceBackedDataSplitter,
    SemiSupervisedDataSplitter,
)
from ._samplers import BatchDistributedSampler
from ._semi_dataloader import SemiSupervisedDataLoader

__all__ = [
    "AnnDataLoader",
    "AnnTorchDataset",
    "CollectionAdapter",
    "ConcatDataLoader",
    "DeviceBackedDataSplitter",
    "SemiSupervisedDataLoader",
    "DataSplitter",
    "SemiSupervisedDataSplitter",
    "BatchDistributedSampler",
    "MappedCollectionDataModule",
    "TileDBDataModule",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
