# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from ml_collections.config_dict import FrozenConfigDict


class attrdict(FrozenConfigDict):
    """A dictionary that allows for attribute-style access.

    Dummy class that allows for backwards compatibility with the previous custom ``attrdict``.
    Inherits from :class:`~ml_collections.config_dict.FrozenConfigDict`.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
