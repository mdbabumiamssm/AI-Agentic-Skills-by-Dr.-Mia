# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

# Simple CoT Demo script (Mock)

def build_cot_prompt(question, method="zero-shot"):
    if method == "zero-shot":
        return f"Question: {question}\nAnswer: Let's think step by step."
    
    elif method == "few-shot":
        examples = """
Q: Roger has 5 tennis balls. He buys 2 more cans of tennis balls. Each can has 3 tennis balls. How many tennis balls does he have now?
A: Roger started with 5 balls. 2 cans of 3 balls each is 6 balls. 5 + 6 = 11. The answer is 11.

Q: The cafeteria had 23 apples. If they used 20 to make lunch and bought 6 more, how many apples do they have?
A: They started with 23. They used 20, so 23 - 20 = 3. They bought 6 more, so 3 + 6 = 9. The answer is 9.
"""
        return f"{examples}\nQ: {question}\nA:"

if __name__ == "__main__":
    q = "A juggler can juggle 16 balls. Half of the balls are golf balls, and half of the golf balls are blue. How many blue golf balls are there?"
    
    print("--- Zero Shot Prompt ---")
    print(build_cot_prompt(q, "zero-shot"))
    print("\n--- Few Shot Prompt ---")
    print(build_cot_prompt(q, "few-shot"))

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
