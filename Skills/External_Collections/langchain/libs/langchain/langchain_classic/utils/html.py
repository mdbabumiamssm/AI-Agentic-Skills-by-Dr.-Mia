# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from langchain_core.utils.html import (
    DEFAULT_LINK_REGEX,
    PREFIXES_TO_IGNORE,
    PREFIXES_TO_IGNORE_REGEX,
    SUFFIXES_TO_IGNORE,
    SUFFIXES_TO_IGNORE_REGEX,
    extract_sub_links,
    find_all_links,
)

__all__ = [
    "DEFAULT_LINK_REGEX",
    "PREFIXES_TO_IGNORE",
    "PREFIXES_TO_IGNORE_REGEX",
    "SUFFIXES_TO_IGNORE",
    "SUFFIXES_TO_IGNORE_REGEX",
    "extract_sub_links",
    "find_all_links",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
