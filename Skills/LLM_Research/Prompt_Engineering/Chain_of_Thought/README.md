<!--
# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

-->

# Chain of Thought (CoT) Prompting

**Category:** LLM Research / Prompt Engineering
**Difficulty:** Beginner

## Overview
Chain of Thought (CoT) prompting improves the reasoning abilities of Large Language Models (LLMs) by encouraging them to generate intermediate reasoning steps before answering a question. 

First proposed in [Wei et al. (2022)](https://arxiv.org/abs/2201.11903), it significantly boosts performance on arithmetic, commonsense, and symbolic reasoning tasks.

## Mechanism
Standard prompting asks for an answer directly:
> Q: Roger has 5 tennis balls. He buys 2 more cans of tennis balls. Each can has 3 tennis balls. How many tennis balls does he have now?
> A: 11.

CoT prompting asks for the reasoning:
> Q: Roger has 5 tennis balls. He buys 2 more cans of tennis balls. Each can has 3 tennis balls. How many tennis balls does he have now?
> A: Roger started with 5 balls. 2 cans of 3 balls each is 6 balls. 5 + 6 = 11. The answer is 11.

## Techniques

### 1. Zero-Shot CoT
Simply appending "Let's think step by step" to the prompt often triggers CoT reasoning without specific examples.
*   **Prompt:** "Solve this problem. Let's think step by step."

### 2. Few-Shot CoT
Providing examples (shots) where the reasoning is explicitly written out.

## Templates

### Arithmetic Template
```text
Q: If I have 3 apples and buy 4 more, how many do I have?
A: I started with 3. I bought 4. 3 + 4 = 7. The answer is 7.

Q: If John has 10 candies and gives 2 to Mary and 3 to Bob, how many does he have left?
A: John starts with 10. He gives 2 to Mary, so 10 - 2 = 8. He gives 3 to Bob, so 8 - 3 = 5. The answer is 5.

Q: {user_question}
A:
```

### Logical Reasoning Template
```text
Q: All men are mortal. Socrates is a man. Is Socrates mortal?
A: Since Socrates is a man, and all men are mortal, it follows that Socrates must be mortal. The answer is Yes.

Q: {user_question}
A:
```

## Best Practices
1.  **Be specific:** Clearly distinguish the "Reasoning" phase from the "Answer" phase.
2.  **Use delimiters:** Use `Reflect:` or `Thought:` to guide the model.
3.  **Validate:** CoT can sometimes produce "hallucinated reasoning" where the steps are wrong but the answer is right (or vice versa).


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->