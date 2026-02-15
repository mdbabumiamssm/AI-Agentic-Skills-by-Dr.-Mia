# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Summarization checker chain for verifying accuracy of text generation.

Chain that tries to verify the accuracy of text generation by splitting it into a
list of facts, then checking if those facts are true or not, and rewriting
the text to make it more truthful. It will repeat this loop until it hits `max_tries` or
gets to a "true" output.
"""

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
