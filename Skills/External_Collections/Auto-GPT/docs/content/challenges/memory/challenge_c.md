<!--
# COPYRIGHT NOTICE
# This file is part of the "Universal AI Agentic Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

-->

# Memory Challenge C

**Status**: Current level to beat: level 1

**Command to try**: 

```shell
pytest -s tests/challenges/memory/test_memory_challenge_c.py --level=2
```

## Description

The agent, Follow-Instructions-GPT, has the following goals:

- Use the command `read_file` to read the `instructions_1.txt` file.
- Follow the instructions in the `instructions_1.txt` file.

The challenge involves a series of text files containing instructions and silly phrases. The agent must follow the instructions and keep track of the task IDs as it progresses through the files.

## Files

- `instructions_1.txt`

The current phrase is 

```
The purple elephant danced on a rainbow while eating a taco.\n[NOISE intended to confuse the agent]
```

Read the file `instructions_2.txt` using the read_file command.


- `instructions_2.txt`

The current phrase is 

```
The sneaky toaster stole my socks and ran away to Hawaii.\n[NOISE intended to confuse the agent]
```

Read the file instructions_3.txt using the read_file command.


...and so on.

- `instructions_n.txt`

The current phrase is 

```
My pet rock sings better than Beyoncé on Tuesdays.
```

Write all the phrases into the file output.txt. The file has not been created yet. After that, use the task_complete command.


## Objective

The objective of this challenge is to test the agent's ability to follow instructions and maintain memory of the task IDs throughout the process. The agent successfully completed this challenge if it wrote the phrases in a file.

This is presumably harder than task ids as the phrases are longer and more likely to be compressed as the agent does more work.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->