# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from agno.agent import Agent
from dotenv import load_dotenv
from typing import List

load_dotenv()


def user_source_selection_run(
    agent: Agent,
    selected_sources: List[int],
) -> str:
    """
    User Source Selection that takes the selected sources indices as input and updates the final confirmed sources.
    Args:
        agent: The agent instance
        selected_sources: The selected sources indices
    Returns:
        Response status
    """
    from services.internal_session_service import SessionService
    session_id = agent.session_id
    session = SessionService.get_session(session_id)
    session_state = session["state"]
    for i, src in enumerate(session_state["search_results"]):
        if (i+1) in selected_sources:
            src["confirmed"] = True
        else:
            src["confirmed"] = False
    SessionService.save_session(session_id, session_state)
    return f"Updated selected sources to {selected_sources}."

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
