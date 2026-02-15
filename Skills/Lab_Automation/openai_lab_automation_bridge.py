# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""
OpenAI Health Stack – Lab Automation Bridge
------------------------------------------
Wraps ExperimentDesigner outputs in a JSON schema compatible with
OpenAI's workflow automation APIs and the Thermo Fisher integration narrative.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any, Dict

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.append(str(ROOT))

from Lab_Automation.Experiment_Design.designer import ExperimentDesigner  # pylint: disable=wrong-import-position


class LabAutomationBridge:
    def __init__(self) -> None:
        self.designer = ExperimentDesigner()

    def synthesize_protocol(self, intent: str, parameters: Dict[str, Any]) -> Dict[str, Any]:
        base_protocol = self.designer.design_experiment(intent, parameters)
        enriched = {
            "intent": intent,
            "schema_version": "openai.lab.2026-01",
            "protocol": base_protocol,
            "handoff": {
                "format": "json",
                "checksum": self._compute_checksum(base_protocol),
                "publish_topic": "openai.health.lab_protocol",
            },
        }
        self._validate(enriched)
        return enriched

    def _compute_checksum(self, protocol: Dict[str, Any]) -> str:
        # Simple checksum for demo (hash of serialized content)
        return hex(abs(hash(json.dumps(protocol, sort_keys=True))))[2:]

    def _validate(self, payload: Dict[str, Any]) -> None:
        if "protocol" not in payload or "handoff" not in payload:
            raise ValueError("Invalid OpenAI lab payload.")
        if "steps" not in payload["protocol"]:
            raise ValueError("Protocol missing steps.")
        if not isinstance(payload["protocol"]["steps"], list):
            raise ValueError("Protocol steps must be a list.")


def _demo() -> None:
    bridge = LabAutomationBridge()
    parameters = {"robot": "Opentrons_OT2", "start_conc": 0, "end_conc": 50, "wells": 8}
    payload = bridge.synthesize_protocol("Dose response gradient", parameters)
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    _demo()

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
