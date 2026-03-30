"""Executable stub for wiring OpenAI Swarm into repo tools."""

import argparse
import json
from pathlib import Path

try:
    from swarm import Agent, Swarm
except ImportError as exc:  # pragma: no cover
    raise SystemExit("Install the `swarm` package before running this demo") from exc

ROOT = Path(__file__).resolve().parents[3]  # project root
TOOLS_REGISTRY = ROOT / "platform" / "adapters" / "runtime_adapter.py"


def load_prompt(path: Path) -> str:
    return path.read_text().strip()


def build_agents():
    triage = Agent(
        name="ClinicalTriage",
        instructions=load_prompt(Path(__file__).with_name("triage_prompt.txt")),
    )
    safety = Agent(
        name="SafetyOfficer",
        instructions=load_prompt(Path(__file__).with_name("safety_prompt.txt")),
    )

    def escalate_if_needed(state):
        if "red_flag" in state.get("summary", ""):
            return safety
        return None

    triage.functions = [escalate_if_needed]
    return triage


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("question", help="Clinical or research question to route through Swarm")
    args = parser.parse_args()

    agent = build_agents()
    swarm = Swarm()
    result = swarm.run(agent=agent, messages=[{"role": "user", "content": args.question}])
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
