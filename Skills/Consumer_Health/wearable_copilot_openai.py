"""
OpenAI Health Stack – Care Copilot
---------------------------------
Transforms wearable and patient-reported data into a schema-enforced
coaching plan similar to ChatGPT Health. 
"""

from __future__ import annotations

import json
from dataclasses import dataclass, asdict
from datetime import datetime
from statistics import mean
from typing import Any, Dict, List, Optional


class SchemaValidationError(Exception):
    """Raised when a generated payload violates the declared schema."""


@dataclass
class CarePlan:
    user_id: str
    generated_at: str
    vitals_summary: Dict[str, Any]
    risk_flags: List[Dict[str, Any]]
    recommendations: List[Dict[str, Any]]

    def to_json(self) -> str:
        return json.dumps(asdict(self), indent=2)


class OpenAIHealthCareCopilot:
    """
    Harmonises wearable telemetry and emits OpenAI-style JSON payloads.
    The schema mirrors public ChatGPT Health demos (summary + actions).
    """

    REQUIRED_VITAL_KEYS = {"resting_heart_rate", "sleep_score", "hrr_variability"}

    def __init__(self) -> None:
        self.stream_buffer: List[Dict[str, Any]] = []

    def ingest_wearable_stream(self, samples: List[Dict[str, Any]]) -> None:
        """Accepts streamed wearable points and keeps the last 30 entries."""
        self.stream_buffer.extend(samples)
        self.stream_buffer = sorted(self.stream_buffer, key=lambda s: s["timestamp"])[-30:]

    def generate_care_plan(self, patient_profile: Dict[str, Any]) -> CarePlan:
        """
        Generates a schema-compliant care plan.
        `patient_profile` contains demographics + optional clinician goals.
        """
        vitals = self._summarise_vitals()
        flags = self._detect_risks(vitals, patient_profile)
        recommendations = self._build_recommendations(flags, patient_profile)

        plan = CarePlan(
            user_id=patient_profile.get("user_id", "anonymous"),
            generated_at=datetime.utcnow().isoformat(),
            vitals_summary=vitals,
            risk_flags=flags,
            recommendations=recommendations,
        )
        self._validate_schema(asdict(plan))
        return plan

    # ------------------
    # Internal helpers
    # ------------------

    def _summarise_vitals(self) -> Dict[str, Any]:
        if not self.stream_buffer:
            return {}

        field_history: Dict[str, List[float]] = {}
        for entry in self.stream_buffer:
            for key, value in entry.items():
                if key == "timestamp":
                    continue
                field_history.setdefault(key, []).append(float(value))

        summary = {}
        for key, values in field_history.items():
            summary[key] = {
                "latest": values[-1],
                "avg_7d": mean(values[-7:]) if len(values) >= 7 else mean(values),
                "trend": values[-1] - values[0],
            }
        return summary

    def _detect_risks(
        self, vitals: Dict[str, Any], patient_profile: Dict[str, Any]
    ) -> List[Dict[str, Any]]:
        risk_flags: List[Dict[str, Any]] = []
        thresholds = patient_profile.get(
            "thresholds",
            {"resting_heart_rate": 5, "sleep_score": -10, "hrr_variability": -15},
        )

        for metric, stats in vitals.items():
            change = stats.get("trend", 0)
            if metric == "resting_heart_rate" and change > thresholds["resting_heart_rate"]:
                risk_flags.append(
                    {
                        "metric": metric,
                        "severity": "high" if change > thresholds["resting_heart_rate"] * 2 else "medium",
                        "message": f"Resting HR increased by {change:.1f} bpm.",
                    }
                )
            if metric == "sleep_score" and change < thresholds["sleep_score"]:
                risk_flags.append(
                    {
                        "metric": metric,
                        "severity": "medium",
                        "message": "Sleep score trending downward more than expected.",
                    }
                )
            if metric == "hrr_variability" and change < thresholds["hrr_variability"]:
                risk_flags.append(
                    {
                        "metric": metric,
                        "severity": "high",
                        "message": "Heart rate variability declining—possible recovery issues.",
                    }
                )
        return risk_flags

    def _build_recommendations(
        self, flags: List[Dict[str, Any]], patient_profile: Dict[str, Any]
    ) -> List[Dict[str, Any]]:
        if not flags:
            return [
                {
                    "goal": "maintain",
                    "action": "Continue current routines; no anomalies detected.",
                    "category": "general",
                }
            ]

        recs: List[Dict[str, Any]] = []
        for flag in flags:
            if flag["metric"] == "resting_heart_rate":
                recs.append(
                    {
                        "goal": "reduce stress load",
                        "action": "Schedule low-intensity recovery day, ensure hydration, and revisit stimulant intake.",
                        "category": "cardio",
                    }
                )
            elif flag["metric"] == "sleep_score":
                recs.append(
                    {
                        "goal": "improve sleep hygiene",
                        "action": "Enforce a fixed bedtime, limit evening screen time, and track caffeine cut-off.",
                        "category": "sleep",
                    }
                )
            else:
                recs.append(
                    {
                        "goal": "boost recovery",
                        "action": "Introduce guided breathing or HRV biofeedback session before bed.",
                        "category": "recovery",
                    }
                )
        return recs

    def _validate_schema(self, payload: Dict[str, Any]) -> None:
        missing = [k for k in CarePlan.__annotations__.keys() if k not in payload]
        if missing:
            raise SchemaValidationError(f"CarePlan payload missing keys: {missing}")
        if not isinstance(payload["risk_flags"], list):
            raise SchemaValidationError("`risk_flags` must be a list.")
        if not isinstance(payload["recommendations"], list):
            raise SchemaValidationError("`recommendations` must be a list.")


def _demo() -> None:
    stream = [
        {"timestamp": datetime(2026, 1, i + 1).isoformat(), "resting_heart_rate": 60 + i, "sleep_score": 85 - i}
        for i in range(7)
    ]
    for entry in stream:
        entry["hrr_variability"] = 80 - (i := stream.index(entry)) * 2

    copilot = OpenAIHealthCareCopilot()
    copilot.ingest_wearable_stream(stream)
    profile = {"user_id": "patient_001", "thresholds": {"resting_heart_rate": 4, "sleep_score": -8, "hrr_variability": -10}}
    plan = copilot.generate_care_plan(profile)
    print(plan.to_json())


if __name__ == "__main__":
    _demo()
