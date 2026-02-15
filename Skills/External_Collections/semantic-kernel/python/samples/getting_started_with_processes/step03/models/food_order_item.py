# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

# Copyright (c) Microsoft. All rights reserved.

from enum import Enum


class FoodItem(str, Enum):
    POTATO_FRIES = "Potato Fries"
    FRIED_FISH = "Fried Fish"
    FISH_SANDWICH = "Fish Sandwich"
    FISH_AND_CHIPS = "Fish & Chips"

    def to_friendly_string(self) -> str:
        return self.value

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
