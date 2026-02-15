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

# Execution Methods

Demonstrates the three execution methods available in the OpenAI Agents SDK: sync, async, and streaming.

## 🎯 What This Demonstrates

- **Runner.run()**: Asynchronous execution for non-blocking operations
- **Runner.run_sync()**: Synchronous execution for simple blocking calls
- **Runner.run_streamed()**: Streaming execution for real-time responses
- **Performance Comparison**: When to use each method

## 🚀 Quick Start

1. **Install OpenAI Agents SDK**:
   ```bash
   pip install openai-agents
   ```

2. **Set up environment**:
   ```bash
   cp ../env.example .env
   # Edit .env and add your OpenAI API key
   ```

3. **Run the agent**:
   ```python
   from agents import Runner
   from agent import root_agent
   
   # Test sync execution
   result = root_agent.sync_execution_example()
   print(result)
   ```

## 💡 Key Concepts

- **Sync Execution**: Blocks until completion, simple to use
- **Async Execution**: Non-blocking, enables concurrency
- **Streaming Execution**: Real-time response processing
- **Use Case Selection**: Choose based on application needs

## 🔗 Next Steps

- [Conversation Management](../4_2_conversation_management/README.md) - Threading and sessions
- [Run Configuration](../4_3_run_configuration/README.md) - Advanced settings


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->