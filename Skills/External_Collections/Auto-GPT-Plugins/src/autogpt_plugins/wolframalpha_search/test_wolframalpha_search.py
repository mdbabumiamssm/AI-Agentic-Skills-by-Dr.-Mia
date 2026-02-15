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
import unittest

import requests

from . import AutoGPTWolframAlphaSearch


class TestAutoGPTWolframAlphaSearch(unittest.TestCase):
    def setUp(self):
        os.environ["WOLFRAMALPHA_APPID"] = "test_appid"
        self.plugin = AutoGPTWolframAlphaSearch()

    def tearDown(self):
        os.environ.pop("WOLFRAMALPHA_APPID", None)

    def test_wolframalpha_search(self):
        query = "2+2"
        try:
            from .wolframalpha_search import _wolframalpha_search
            _wolframalpha_search(query)
        except requests.exceptions.HTTPError as e:
            self.assertEqual(e.response.status_code, 401)


if __name__ == "__main__":
    unittest.main()

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
