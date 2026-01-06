"""
Autonomous Lab Controller (2026 Skills)

This agent acts as a 'Compiler' for Biology.
Input: High-level natural language intent.
Output: Low-level JSON protocol for robotic liquid handlers (e.g., Opentrons/Tecan).

Features:
- Parameter extraction (Volume, Reagents, Wells).
- Safety checks (Volume limits).
- Protocol serialization.
"""

import json
import re

class RobotProtocolCompiler:
    def __init__(self, robot_model="Opentrons_OT2"):
        self.robot_model = robot_model
        self.max_volume_ul = 1000  # P1000 pipette
        self.deck_layout = {
            "1": "tiprack_1000ul",
            "2": "96_wellplate_flat",
            "3": "reagent_reservoir"
        }

    def parse_intent(self, user_instruction: str) -> dict:
        """
        Extracts semantic actions from natural language.
        In production, this would use an LLM. Here, we use Regex for the demo.
        """
        protocol = {
            "meta": {
                "robot": self.robot_model,
                "task": user_instruction
            },
            "steps": []
        }

        # Regex patterns for common lab actions
        # "Transfer 100ul of Water to A1"
        transfer_pattern = r"Transfer (\d+)ul of (\w+) to ([A-H]\d+)"
        
        matches = re.findall(transfer_pattern, user_instruction, re.IGNORECASE)
        
        for vol, reagent, dest in matches:
            step = {
                "action": "transfer",
                "volume": int(vol),
                "source": reagent,
                "destination": dest,
                "tool": "p1000" if int(vol) > 20 else "p20"
            }
            protocol["steps"].append(step)

        # "Mix A1 5 times"
        mix_pattern = r"Mix ([A-H]\d+) (\d+) times"
        mix_matches = re.findall(mix_pattern, user_instruction, re.IGNORECASE)
        
        for well, times in mix_matches:
            step = {
                "action": "mix",
                "well": well,
                "repetitions": int(times),
                "volume": 50 # Default mix volume
            }
            protocol["steps"].append(step)
            
        return protocol

    def validate_protocol(self, protocol_data: dict) -> list:
        """
        Safety checks for the robot.
        """
        errors = []
        for i, step in enumerate(protocol_data["steps"]):
            if step["action"] == "transfer":
                if step["volume"] > self.max_volume_ul:
                    errors.append(f"Step {i+1}: Volume {step['volume']}ul exceeds pipette limit ({self.max_volume_ul}ul).")
        return errors

    def generate_json(self, user_instruction: str):
        print(f"--- Compiling Protocol: '{user_instruction}' ---")
        
        # 1. Parse
        data = self.parse_intent(user_instruction)
        
        # 2. Validate
        errors = self.validate_protocol(data)
        if errors:
            print("COMPILATION FAILED:")
            for e in errors:
                print(f"  - {e}")
            return None
            
        # 3. Output
        json_output = json.dumps(data, indent=2)
        print("Protocol Generated Successfully:")
        print(json_output)
        return json_output

if __name__ == "__main__":
    bot = RobotProtocolCompiler()
    
    # Test 1: Valid
    instruction = "Transfer 200ul of Buffer to A1. Then transfer 50ul of Sample to A1. Finally Mix A1 10 times."
    bot.generate_json(instruction)
    
    print("\n" + "="*40 + "\n")
    
    # Test 2: Error (Volume too high)
    bot.generate_json("Transfer 1500ul of Media to B2.")
