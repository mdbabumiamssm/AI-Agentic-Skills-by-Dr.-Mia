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

from promptflow import PFClient

pf_client = PFClient()

inputs = {
    "text": "What would you have left if you spent $3 when you only had $2 to begin with"
}  # The inputs of the flow.
flow_result = pf_client.test(flow="perform_math", inputs=inputs)
print(f"Flow outputs: {flow_result}")

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
