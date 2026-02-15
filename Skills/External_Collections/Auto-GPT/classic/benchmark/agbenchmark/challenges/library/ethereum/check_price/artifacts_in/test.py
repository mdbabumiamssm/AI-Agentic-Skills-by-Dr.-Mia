# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import re

from .sample_code import get_ethereum_price


def test_get_ethereum_price() -> None:
    # Read the Ethereum price from the file
    with open("eth_price.txt", "r") as file:
        eth_price = file.read().strip()

    # Validate that the eth price is all digits
    pattern = r"^\d+$"
    matches = re.match(pattern, eth_price) is not None
    assert (
        matches
    ), f"AssertionError: Ethereum price should be all digits, but got {eth_price}"

    # Get the current price of Ethereum
    real_eth_price = get_ethereum_price()

    # Convert the eth price to a numerical value for comparison
    eth_price_value = float(eth_price)
    real_eth_price_value = float(real_eth_price)

    # Check if the eth price is within $50 of the actual Ethereum price
    assert abs(real_eth_price_value - eth_price_value) <= 50, (
        "AssertionError: Ethereum price is not within $50 of the actual Ethereum price "
        f"(Provided price: ${eth_price}, Real price: ${real_eth_price})"
    )

    print("Matches")


if __name__ == "__main__":
    test_get_ethereum_price()

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
