# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import unittest

from .url_shortener import retrieve_url, shorten_url


class TestURLShortener(unittest.TestCase):
    def test_url_retrieval(self):
        # Shorten the URL to get its shortened form
        shortened_url = shorten_url("https://www.example.com")

        # Retrieve the original URL using the shortened URL directly
        retrieved_url = retrieve_url(shortened_url)

        self.assertEqual(
            retrieved_url,
            "https://www.example.com",
            "Retrieved URL does not match the original!",
        )


if __name__ == "__main__":
    unittest.main()

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
