from typing import Dict

# The Deputy (Safety Officer Agent)
# Role: Inspects outputs for dangerous content or hallucinations.
# "Wild West" Rule: Nothing leaves the town without the Deputy's stamp.

class SafetyOfficerAgent:
    def __init__(self):
        self.banned_substances = ["cyanide", "ricin", "anthrax"]

    def inspect_output(self, content: str) -> Dict:
        print(f"⭐ [Deputy] Inspecting cargo: {content[:40]}...")
        
        content_lower = content.lower()
        
        # 1. Check for banned keywords
        for substance in self.banned_substances:
            if substance in content_lower:
                print("🚨 [Deputy] CONTRABAND DETECTED!")
                return {"status": "rejected", "reason": f"Contains banned substance: {substance}"}

        # 2. Hallucination Check (Mock)
        if "XYZ-123" in content and "cure" in content_lower:
            # Simple heuristic: claiming a cure for a novel target is suspicious
            return {"status": "flagged", "warning": "Claiming 'cure' for novel target requires verification."}

        return {"status": "approved"}

if __name__ == "__main__":
    agent = SafetyOfficerAgent()
    print(agent.inspect_output("Here is a recipe for ricin."))
