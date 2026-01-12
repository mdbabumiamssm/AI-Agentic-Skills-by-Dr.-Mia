"""
Anthropic Health Stack – Pharmacovigilance Monitor
-------------------------------------------------
Analyzes adverse event feeds and prioritizes them following Anthropic
\"coworker\" conventions (trace text + structured JSON event).
"""

from __future__ import annotations

import json
from dataclasses import dataclass, asdict
from datetime import datetime
from typing import Dict, List


@dataclass
class SafetySignal:
    case_id: str
    product: str
    event_term: str
    seriousness: str
    narrative: str
    recommended_action: str
    generated_at: str


class PharmacovigilanceMonitor:
    SERIOUS_TERMS = {"anaphylaxis", "death", "hospitalisation"}

    def triage(self, cases: List[Dict[str, str]]) -> List[Dict[str, str]]:
        signals: List[Dict[str, str]] = []
        for case in cases:
            seriousness = self._score_seriousness(case["event_term"])
            action = "Immediate medical review" if seriousness == "high" else "Add to monitoring queue"
            entry = SafetySignal(
                case_id=case["case_id"],
                product=case["product"],
                event_term=case["event_term"],
                seriousness=seriousness,
                narrative=self._build_trace(case, seriousness),
                recommended_action=action,
                generated_at=datetime.utcnow().isoformat(),
            )
            signals.append(asdict(entry))
        return signals

    def _score_seriousness(self, event_term: str) -> str:
        lowered = event_term.lower()
        return "high" if any(term in lowered for term in self.SERIOUS_TERMS) else "moderate"

    def _build_trace(self, case: Dict[str, str], seriousness: str) -> str:
        return (
            f"<thinking>Evaluating case {case['case_id']} ({case['product']}).</thinking>\n"
            f"<analysis>Event '{case['event_term']}' mapped to seriousness '{seriousness}'.</analysis>\n"
            f"<decision>Recommend { 'urgent escalation' if seriousness == 'high' else 'routine monitoring'}.</decision>"
        )


def _demo() -> None:
    monitor = PharmacovigilanceMonitor()
    payload = monitor.triage(
        [
            {"case_id": "AE-1", "product": "RespiroX", "event_term": "anaphylaxis"},
            {"case_id": "AE-2", "product": "RespiroX", "event_term": "mild rash"},
        ]
    )
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    _demo()
