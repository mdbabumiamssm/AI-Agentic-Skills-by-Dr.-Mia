"""
Anthropic Health Stack – Regulatory Coworker
-------------------------------------------
Drafts responses to regulatory agencies with structured justifications.
"""

from __future__ import annotations

import json
from datetime import datetime
from typing import Dict, List


class RegulatoryCoworker:
    def __init__(self) -> None:
        self.ctd_sections = {"2.3.S", "2.4", "2.5"}

    def prepare_response(self, query: Dict[str, str], evidence: List[str]) -> Dict[str, str]:
        section = query.get("section", "2.3.S")
        if section not in self.ctd_sections:
            raise ValueError("Unsupported CTD section.")

        citations = [f"doc://{i}" for i, _ in enumerate(evidence)]

        narrative = (
            f"<thinking>Clarifying {section} question: {query['question']}</thinking>\n"
            f"<analysis>Reviewed {len(evidence)} evidence files.</analysis>\n"
            "<decision>Drafting response referencing cited documents.</decision>"
        )

        draft = {
            "section": section,
            "query_id": query.get("id", "unknown"),
            "response": f"Per section {section}, the manufacturing process ...",
            "citations": citations,
            "trace": narrative,
            "generated_at": datetime.utcnow().isoformat(),
        }
        return draft


def _demo() -> None:
    coworker = RegulatoryCoworker()
    query = {"id": "FDA-REQ-9", "section": "2.3.S", "question": "Clarify impurity profile."}
    payload = coworker.prepare_response(query, ["/tmp/impurity_report.pdf"])
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    _demo()
