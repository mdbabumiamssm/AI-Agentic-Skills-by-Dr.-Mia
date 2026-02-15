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

from semantic_kernel.processes.kernel_process.kernel_process_edge import KernelProcessEdge
from semantic_kernel.processes.kernel_process.kernel_process_function_target import KernelProcessFunctionTarget
from semantic_kernel.processes.process_message_factory import ProcessMessageFactory


def test_create_from_edge():
    """Test initialization of KernelProcessEdge with valid input."""
    source_step_id = "step_1"
    output_target = KernelProcessFunctionTarget(
        step_id="step_2", function_name="process_data", parameter_name="input_data", target_event_id="event_1"
    )
    edge = KernelProcessEdge(source_step_id=source_step_id, output_target=output_target)

    process_message = ProcessMessageFactory.create_from_edge(edge, "data")

    assert process_message.source_id == source_step_id

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
