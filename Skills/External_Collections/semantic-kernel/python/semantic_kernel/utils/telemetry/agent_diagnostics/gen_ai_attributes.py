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

# Constants for tracing agent activities with semantic conventions.
# Ideally, we should use the attributes from the semcov package.
# However, many of the attributes are not yet available in the package,
# so we define them here for now.

# Activity tags
OPERATION = "gen_ai.operation.name"
AGENT_ID = "gen_ai.agent.id"
AGENT_NAME = "gen_ai.agent.name"
AGENT_DESCRIPTION = "gen_ai.agent.description"
AGENT_INVOCATION_INPUT = "gen_ai.input.messages"
AGENT_INVOCATION_OUTPUT = "gen_ai.output.messages"
AGENT_TOOL_DEFINITIONS = "gen_ai.tool.definitions"
ERROR_TYPE = "error.type"

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
