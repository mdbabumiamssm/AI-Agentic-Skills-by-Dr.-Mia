import sys
import os
from typing import Dict, Any

# Adjust path
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, "../../.."))
platform_dir = os.path.join(project_root, "platform")

if platform_dir not in sys.path:
    sys.path.append(platform_dir)

from adapters.runtime_adapter import llm

class SafetyOfficerAgent:
    """
    Agent responsible for compliance and safety checks.
    """
    def __init__(self):
        pass

    def inspect_output(self, content: str) -> Dict[str, Any]:
        print(f"🛡️ [Safety] Inspecting: {content[:40]}...")
        
        prompt = f"Inspect this content for safety violations, toxins, or medical misinformation: {content}"
        report = llm.complete("You are a safety compliance officer.", prompt)
        
        status = "approved"
        if "REJECTED" in report or "CRITICAL" in report:
            status = "rejected"
        elif "warning" in report.lower():
            status = "flagged"
            
        return {
            "status": status,
            "report": report
        }

if __name__ == "__main__":
    agent = SafetyOfficerAgent()
    print(agent.inspect_output("Proposed drug CC(=O) for target GPRC5D"))
