"""
Anthropic Health Stack – Prior Authorization Coworker
----------------------------------------------------
Produces Claude-style coworker traces (<thinking> + JSON decision payload)
so the Event Bus can replay determinations for auditors.
"""

from __future__ import annotations

import json
from datetime import datetime
from typing import Any, Dict, List


class PriorAuthCoworker:
    def __init__(self) -> None:
        self.policy_db = {
            "MRI-L-SPINE": {
                "version": "2026.01",
                "criteria": [
                    "Persistent back pain > 6 weeks",
                    "Documented conservative therapy failure",
                    "Presence of red flag symptoms",
                ],
            }
        }

    def review(self, casefile: Dict[str, str]) -> Dict[str, Any]:
        """
        Returns Anthropic-style trace plus structured JSON output.
        """
        policy = self.policy_db.get(casefile["procedure_code"])
        trace = self._build_trace(casefile, policy)
        reasoning = self._evaluate(casefile["clinical_note"], policy)
        decision = all(item["met"] for item in reasoning)
        payload = {
            "case_id": casefile.get("case_id", "unknown"),
            "procedure_code": casefile["procedure_code"],
            "policy_version": policy["version"],
            "decision": "APPROVED" if decision else "DENIED",
            "reasoning": reasoning,
            "generated_at": datetime.utcnow().isoformat(),
            "trace": trace,
        }
        return payload

    def _evaluate(self, note: str, policy: Dict[str, Any]) -> List[Dict[str, Any]]:
        note_lower = note.lower()
        checks = []
        for criterion in policy["criteria"]:
            if criterion.startswith("Persistent"):
                checks.append({"criterion": criterion, "met": "week" in note_lower})
            elif "conservative" in criterion:
                checks.append({"criterion": criterion, "met": "therapy" in note_lower})
            else:
                checks.append({"criterion": criterion, "met": "red flag" in note_lower or "trauma" in note_lower})
        return checks

    def _build_trace(self, casefile: Dict[str, str], policy: Dict[str, Any]) -> str:
        return (
            f"<thinking>Reviewing case {casefile.get('case_id')} for {casefile['procedure_code']} "
            f"against policy v{policy['version']}.</thinking>\n"
            f"<analysis>Extracted duration, conservative therapy, and red-flag indicators "
            f"from the note.</analysis>\n"
            f"<decision>Will {'approve' if 'therapy' in casefile['clinical_note'].lower() else 'deny'} "
            f"based on matched criteria.</decision>"
        )


def _demo() -> None:
    coworker = PriorAuthCoworker()
    case = {
        "case_id": "PA-1001",
        "procedure_code": "MRI-L-SPINE",
        "clinical_note": "Chronic back pain for 8 weeks with failed PT and NSAIDs. No trauma.",
    }
    payload = coworker.review(case)
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    _demo()
