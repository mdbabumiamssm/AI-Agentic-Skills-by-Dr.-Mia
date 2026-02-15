# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from langchain_fireworks.chat_models import ChatFireworks
from langchain_fireworks.embeddings import FireworksEmbeddings
from langchain_fireworks.llms import Fireworks
from langchain_fireworks.version import __version__

__all__ = [
    "ChatFireworks",
    "Fireworks",
    "FireworksEmbeddings",
    "__version__",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
