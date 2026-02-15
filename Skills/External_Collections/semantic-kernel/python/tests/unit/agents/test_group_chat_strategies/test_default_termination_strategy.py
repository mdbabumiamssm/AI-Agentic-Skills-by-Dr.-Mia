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


from semantic_kernel.agents.strategies.termination.default_termination_strategy import DefaultTerminationStrategy


async def test_should_agent_terminate_():
    strategy = DefaultTerminationStrategy(maximum_iterations=2)
    result = await strategy.should_agent_terminate(None, [])
    assert not result

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
