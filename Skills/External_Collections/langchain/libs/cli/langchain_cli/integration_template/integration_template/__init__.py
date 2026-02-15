# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from importlib import metadata

from __module_name__.chat_models import Chat__ModuleName__
from __module_name__.document_loaders import __ModuleName__Loader
from __module_name__.embeddings import __ModuleName__Embeddings
from __module_name__.retrievers import __ModuleName__Retriever
from __module_name__.toolkits import __ModuleName__Toolkit
from __module_name__.tools import __ModuleName__Tool
from __module_name__.vectorstores import __ModuleName__VectorStore

try:
    __version__ = metadata.version(__package__)
except metadata.PackageNotFoundError:
    # Case where package metadata is not available.
    __version__ = ""
del metadata  # optional, avoids polluting the results of dir(__package__)

__all__ = [
    "Chat__ModuleName__",
    "__ModuleName__VectorStore",
    "__ModuleName__Embeddings",
    "__ModuleName__Loader",
    "__ModuleName__Retriever",
    "__ModuleName__Toolkit",
    "__ModuleName__Tool",
    "__version__",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
